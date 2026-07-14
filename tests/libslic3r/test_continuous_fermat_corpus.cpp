#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ContinuousFermat.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace Slic3r;

namespace {

namespace stdfs = std::filesystem;
using Clock = std::chrono::steady_clock;

struct CorpusEntry
{
    std::string id;
    stdfs::path path;
    std::string source_url;
    double scale { 1.0 };
};

struct LayerResult
{
    size_t index { 0 };
    double z { 0.0 };
    size_t islands { 0 };
    size_t holes { 0 };
    double area_mm2 { 0.0 };
    double overlap_mm2 { 0.0 };
    double generation_ms { 0.0 };
    double validation_ms { 0.0 };
    double layer_elapsed_ms { 0.0 };
    size_t path_points { 0 };
    size_t reused_from_layer { std::numeric_limits<size_t>::max() };
    bool eligible { false };
    bool attempted { false };
    bool cache_hit { false };
    bool passed { false };
    bool quality_ok { false };
    std::string category;
    std::string reason;
    ContinuousFermat::PathValidation validation;
};

using CanonicalPoint = std::pair<coord_t, coord_t>;
using CanonicalRing = std::vector<CanonicalPoint>;

struct CanonicalGeometry
{
    CanonicalRing outer;
    std::vector<CanonicalRing> holes;

    bool operator==(const CanonicalGeometry &other) const { return outer == other.outer && holes == other.holes; }
};

struct CachedLayerResult
{
    CanonicalGeometry geometry;
    LayerResult result;
    size_t source_layer { 0 };
};

using LayerCache = std::unordered_map<uint64_t, std::vector<CachedLayerResult>>;

const CachedLayerResult *find_cached_layer(const LayerCache &cache, const uint64_t hash, const CanonicalGeometry &geometry)
{
    const auto bucket = cache.find(hash);
    if (bucket == cache.end())
        return nullptr;
    for (const CachedLayerResult &candidate : bucket->second)
        if (candidate.geometry == geometry)
            return &candidate;
    return nullptr;
}

size_t least_rotation(const CanonicalRing &ring)
{
    const size_t count = ring.size();
    if (count < 2)
        return 0;
    size_t first = 0;
    size_t second = 1;
    size_t offset = 0;
    while (first < count && second < count && offset < count) {
        const CanonicalPoint &a = ring[(first + offset) % count];
        const CanonicalPoint &b = ring[(second + offset) % count];
        if (a == b) {
            ++offset;
            continue;
        }
        if (a > b) {
            first += offset + 1;
            if (first <= second)
                first = second + 1;
        } else {
            second += offset + 1;
            if (second <= first)
                second = first + 1;
        }
        offset = 0;
    }
    return std::min(first, second) % count;
}

CanonicalRing rotate_ring(const CanonicalRing &ring, const size_t start)
{
    CanonicalRing result;
    result.reserve(ring.size());
    for (size_t i = 0; i < ring.size(); ++i)
        result.emplace_back(ring[(start + i) % ring.size()]);
    return result;
}

CanonicalRing canonical_ring(const Polygon &polygon)
{
    CanonicalRing forward;
    forward.reserve(polygon.points.size());
    for (const Point &point : polygon.points)
        forward.emplace_back(point.x(), point.y());
    if (forward.empty())
        return forward;
    CanonicalRing reverse(forward.rbegin(), forward.rend());
    forward = rotate_ring(forward, least_rotation(forward));
    reverse = rotate_ring(reverse, least_rotation(reverse));
    return std::min(forward, reverse);
}

CanonicalGeometry canonical_geometry(const ExPolygon &polygon)
{
    CanonicalGeometry result;
    result.outer = canonical_ring(polygon.contour);
    result.holes.reserve(polygon.holes.size());
    for (const Polygon &hole : polygon.holes)
        result.holes.emplace_back(canonical_ring(hole));
    std::sort(result.holes.begin(), result.holes.end());
    return result;
}

uint64_t canonical_hash(const CanonicalGeometry &geometry)
{
    uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&hash](const uint64_t value) {
        hash ^= value;
        hash *= 1099511628211ULL;
    };
    const auto mix_ring = [&mix](const CanonicalRing &ring) {
        mix(uint64_t(ring.size()));
        for (const CanonicalPoint &point : ring) {
            mix(uint64_t(int64_t(point.first)));
            mix(uint64_t(int64_t(point.second)));
        }
    };
    mix_ring(geometry.outer);
    mix(0xa5a5a5a5a5a5a5a5ULL);
    mix(uint64_t(geometry.holes.size()));
    for (const CanonicalRing &hole : geometry.holes) {
        mix(0x5a5a5a5a5a5a5a5aULL);
        mix_ring(hole);
    }
    return hash;
}

double elapsed_ms(const Clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

std::string env_string(const char *name, const std::string &fallback = {})
{
    const char *value = std::getenv(name);
    return value == nullptr || *value == '\0' ? fallback : value;
}

double env_double(const char *name, double fallback)
{
    const std::string value = env_string(name);
    if (value.empty())
        return fallback;
    try {
        return std::stod(value);
    } catch (...) {
        return fallback;
    }
}

size_t env_size(const char *name, size_t fallback)
{
    const std::string value = env_string(name);
    if (value.empty())
        return fallback;
    try {
        return size_t(std::stoull(value));
    } catch (...) {
        return fallback;
    }
}

bool env_bool(const char *name, bool fallback)
{
    std::string value = env_string(name);
    if (value.empty())
        return fallback;
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    return value == "1" || value == "true" || value == "yes" || value == "on";
}

std::string trim_copy(const std::string &value)
{
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) { return std::isspace(c); });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) { return std::isspace(c); }).base();
    return first < last ? std::string(first, last) : std::string {};
}

std::vector<std::string> split(const std::string &value, const char delimiter)
{
    std::vector<std::string> fields;
    size_t begin = 0;
    for (;;) {
        const size_t end = value.find(delimiter, begin);
        fields.emplace_back(value.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos)
            return fields;
        begin = end + 1;
    }
}

std::unordered_set<std::string> parse_model_ids(const std::string &value)
{
    std::unordered_set<std::string> ids;
    if (trim_copy(value).empty())
        return ids;
    for (const std::string &field : split(value, ',')) {
        const std::string id = trim_copy(field);
        if (id.empty())
            throw std::runtime_error("CONTINUOUS_FERMAT_CORPUS_MODEL_IDS contains an empty model ID");
        if (!ids.insert(id).second)
            throw std::runtime_error("CONTINUOUS_FERMAT_CORPUS_MODEL_IDS contains duplicate model ID: " + id);
    }
    return ids;
}

using LayerRange = std::pair<size_t, size_t>;

size_t parse_layer_index(const std::string &value)
{
    if (value.empty() || !std::all_of(value.begin(), value.end(), [](unsigned char c) { return std::isdigit(c); }))
        throw std::runtime_error("invalid zero-based layer index: " + value);
    try {
        const unsigned long long parsed = std::stoull(value);
        if (parsed > std::numeric_limits<size_t>::max())
            throw std::out_of_range("layer index exceeds size_t");
        return size_t(parsed);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid zero-based layer index: " + value);
    }
}

std::vector<LayerRange> parse_layer_ranges(const std::string &value)
{
    std::vector<LayerRange> ranges;
    for (const std::string &raw_field : split(value, ',')) {
        const std::string field = trim_copy(raw_field);
        if (field.empty())
            throw std::runtime_error("layer selector contains an empty index or range");
        const size_t hyphen = field.find('-');
        if (hyphen == std::string::npos) {
            const size_t index = parse_layer_index(field);
            ranges.emplace_back(index, index);
            continue;
        }
        if (field.find('-', hyphen + 1) != std::string::npos)
            throw std::runtime_error("invalid layer range: " + field);
        const size_t first = parse_layer_index(trim_copy(field.substr(0, hyphen)));
        const size_t last = parse_layer_index(trim_copy(field.substr(hyphen + 1)));
        if (last < first)
            throw std::runtime_error("layer range ends before it starts: " + field);
        ranges.emplace_back(first, last);
    }
    std::sort(ranges.begin(), ranges.end());
    std::vector<LayerRange> merged;
    for (const LayerRange &range : ranges) {
        if (merged.empty() || (merged.back().second != std::numeric_limits<size_t>::max() && range.first > merged.back().second + 1))
            merged.emplace_back(range);
        else
            merged.back().second = std::max(merged.back().second, range.second);
    }
    return merged;
}

struct LayerFilter
{
    bool active { false };
    std::vector<LayerRange> global_ranges;
    std::unordered_map<std::string, std::vector<LayerRange>> model_ranges;

    bool selects(const std::string &model_id, const size_t layer_index) const
    {
        if (!active)
            return true;
        const std::vector<LayerRange> *ranges = &global_ranges;
        if (!model_ranges.empty()) {
            const auto found = model_ranges.find(model_id);
            if (found == model_ranges.end())
                return false;
            ranges = &found->second;
        }
        return std::any_of(ranges->begin(), ranges->end(), [layer_index](const LayerRange &range) {
            return layer_index >= range.first && layer_index <= range.second;
        });
    }
};

LayerFilter parse_layer_filter(const std::string &value)
{
    LayerFilter filter;
    const std::string trimmed = trim_copy(value);
    if (trimmed.empty())
        return filter;
    filter.active = true;
    if (trimmed.find('=') == std::string::npos) {
        if (trimmed.find(';') != std::string::npos)
            throw std::runtime_error("global layer selector uses commas, not semicolons");
        filter.global_ranges = parse_layer_ranges(trimmed);
        return filter;
    }
    for (const std::string &raw_clause : split(trimmed, ';')) {
        const std::string clause = trim_copy(raw_clause);
        const size_t equals = clause.find('=');
        if (equals == std::string::npos || clause.find('=', equals + 1) != std::string::npos)
            throw std::runtime_error("per-model layer selector clause must be MODEL_ID=INDICES: " + clause);
        const std::string id = trim_copy(clause.substr(0, equals));
        const std::string ranges = trim_copy(clause.substr(equals + 1));
        if (id.empty() || ranges.empty())
            throw std::runtime_error("per-model layer selector needs a nonempty model ID and layer list");
        if (!filter.model_ranges.emplace(id, parse_layer_ranges(ranges)).second)
            throw std::runtime_error("per-model layer selector contains duplicate model ID: " + id);
    }
    return filter;
}

std::string json_escape(const std::string &value)
{
    std::ostringstream out;
    for (const unsigned char c : value) {
        switch (c) {
        case '\\': out << "\\\\"; break;
        case '"': out << "\\\""; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (c < 0x20)
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
            else
                out << char(c);
        }
    }
    return out.str();
}

std::string csv_escape(const std::string &value)
{
    if (value.find_first_of(",\"\r\n") == std::string::npos)
        return value;
    std::string escaped = value;
    size_t pos = 0;
    while ((pos = escaped.find('"', pos)) != std::string::npos) {
        escaped.insert(pos, 1, '"');
        pos += 2;
    }
    return '"' + escaped + '"';
}

std::vector<std::string> split_tabs(const std::string &line)
{
    std::vector<std::string> fields;
    size_t begin = 0;
    for (;;) {
        const size_t end = line.find('\t', begin);
        fields.emplace_back(line.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos)
            break;
        begin = end + 1;
    }
    return fields;
}

std::vector<CorpusEntry> read_manifest(const stdfs::path &manifest)
{
    std::ifstream input(manifest);
    if (!input)
        throw std::runtime_error("cannot open manifest: " + manifest.string());

    std::vector<CorpusEntry> entries;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty() || line[0] == '#')
            continue;
        const std::vector<std::string> fields = split_tabs(line);
        CorpusEntry entry;
        if (fields.size() == 1) {
            entry.path = stdfs::u8path(fields[0]);
            entry.id = entry.path.stem().u8string() + "-" + std::to_string(line_number);
        } else {
            entry.id = fields[0];
            entry.path = stdfs::u8path(fields[1]);
            if (fields.size() > 2)
                entry.source_url = fields[2];
            if (fields.size() > 3 && !fields[3].empty())
                entry.scale = std::stod(fields[3]);
        }
        if (entry.path.is_relative())
            entry.path = manifest.parent_path() / entry.path;
        entries.emplace_back(std::move(entry));
    }
    return entries;
}

double area_mm2(const ExPolygons &polygons)
{
    return std::abs(area(polygons)) * SCALING_FACTOR * SCALING_FACTOR;
}

double area_mm2(const ExPolygon &polygon)
{
    return std::abs(polygon.area()) * SCALING_FACTOR * SCALING_FACTOR;
}

std::string validation_category(const ContinuousFermat::PathValidation &v)
{
    const auto has = [&v](const char *token) { return v.reason.find(token) != std::string::npos; };
    // Classify from the validator's own reason string rather than duplicating
    // its feature-size-adjusted thresholds here.
    if (has("metadata") || has("multiplier"))
        return "metadata";
    if (has("not exactly closed") || !v.closed)
        return "open_path";
    if (has("crossings="))
        return "crossing";
    if (has("turnbacks="))
        return "turnback";
    if (has("centerline_containment=") || has(" outside="))
        return "outside_or_containment";
    if (has("redeposition=") || has("bead_overlap_pairs="))
        return "redeposition";
    if (has("exact_coverage="))
        return "exact_coverage";
    if (has(" coverage="))
        return "coverage";
    if (has(" material="))
        return "material";
    if (has("width") || has("cross-section") || has("safe range") || has("swept extrusion footprint"))
        return "flow_safety";
    if (has("printable area"))
        return "invalid_area";
    return "validation";
}

void write_layer_json(std::ostream &out, const CorpusEntry &entry, const LayerResult &layer)
{
    const auto &v = layer.validation;
    out << std::setprecision(10)
        << "{\"record\":\"layer\",\"id\":\"" << json_escape(entry.id)
        << "\",\"path\":\"" << json_escape(entry.path.u8string())
        << "\",\"layer\":" << layer.index << ",\"z_mm\":" << layer.z
        << ",\"islands\":" << layer.islands << ",\"holes\":" << layer.holes
        << ",\"area_mm2\":" << layer.area_mm2 << ",\"overlap_mm2\":" << layer.overlap_mm2
        << ",\"eligible\":" << (layer.eligible ? "true" : "false")
        << ",\"attempted\":" << (layer.attempted ? "true" : "false")
        << ",\"passed\":" << (layer.passed ? "true" : "false")
        << ",\"quality_ok\":" << (layer.quality_ok ? "true" : "false")
        << ",\"category\":\"" << json_escape(layer.category) << "\",\"reason\":\"" << json_escape(layer.reason)
        << "\",\"generation_ms\":" << layer.generation_ms << ",\"validation_ms\":" << layer.validation_ms
        << ",\"layer_elapsed_ms\":" << layer.layer_elapsed_ms << ",\"cache_hit\":" << (layer.cache_hit ? "true" : "false")
        << ",\"reused_from_layer\":";
    if (layer.cache_hit)
        out << layer.reused_from_layer;
    else
        out << "null";
    out
        << ",\"path_points\":" << layer.path_points << ",\"closed\":" << (v.closed ? "true" : "false")
        << ",\"emittable\":" << (v.emittable ? "true" : "false")
        << ",\"exact_coverage\":" << v.exact_coverage_ratio << ",\"coverage\":" << v.coverage_ratio
        << ",\"outside\":" << v.outside_ratio << ",\"material\":" << v.material_ratio
        << ",\"redeposition\":" << v.redeposition_ratio << ",\"containment\":" << v.containment_violations
        << ",\"crossings\":" << v.crossings << ",\"spacing\":" << v.spacing_violations
        << ",\"bead_overlaps\":" << v.bead_overlap_violations << ",\"turnbacks\":" << v.turnback_violations
        << ",\"thin_feature_class\":\"" << json_escape(v.thin_feature_class)
        << "\",\"thin_feature_threshold_mm\":" << v.thin_feature_threshold_mm << "}\n";
}

void write_layer_csv(std::ostream &out, const CorpusEntry &entry, const LayerResult &layer)
{
    const auto &v = layer.validation;
    out << std::setprecision(10) << "layer," << csv_escape(entry.id) << ',' << csv_escape(entry.path.u8string()) << ',' << layer.index << ','
        << layer.z << ',' << layer.islands << ',' << layer.holes << ',' << layer.area_mm2 << ',' << layer.overlap_mm2 << ',' << layer.eligible << ','
        << layer.attempted << ',' << layer.passed << ',' << layer.quality_ok << ',' << csv_escape(layer.category) << ',' << csv_escape(layer.reason) << ',' << layer.generation_ms
        << ',' << layer.validation_ms << ',' << layer.path_points << ',' << v.closed << ',' << v.exact_coverage_ratio << ',' << v.coverage_ratio << ','
        << v.outside_ratio << ',' << v.material_ratio << ',' << v.redeposition_ratio << ',' << v.containment_violations << ',' << v.crossings << ','
        << v.spacing_violations << ',' << v.bead_overlap_violations << ',' << v.turnback_violations << ',' << csv_escape(v.thin_feature_class) << ','
        << v.thin_feature_threshold_mm << ",,,,,," << layer.layer_elapsed_ms << ','
        << layer.cache_hit << ',';
    if (layer.cache_hit)
        out << layer.reused_from_layer;
    out << ",,\n";
}

} // namespace

TEST_CASE("Continuous Fermat corpus canonical geometry is exact and hole-aware", "[ContinuousFermat][corpus-canonical]")
{
    const auto point = [](const double x, const double y) { return Point::new_scale(x, y); };
    ExPolygon first({ point(0.0, 0.0), point(12.0, 0.0), point(12.0, 10.0), point(0.0, 10.0) });
    first.holes.emplace_back(Points { point(2.0, 2.0), point(2.0, 4.0), point(4.0, 4.0), point(4.0, 2.0) });
    first.holes.emplace_back(Points { point(7.0, 5.0), point(7.0, 8.0), point(10.0, 8.0), point(10.0, 5.0) });

    ExPolygon reordered = first;
    std::rotate(reordered.contour.points.begin(), reordered.contour.points.begin() + 2, reordered.contour.points.end());
    std::reverse(reordered.contour.points.begin(), reordered.contour.points.end());
    std::swap(reordered.holes[0], reordered.holes[1]);
    std::rotate(reordered.holes[0].points.begin(), reordered.holes[0].points.begin() + 1, reordered.holes[0].points.end());
    std::reverse(reordered.holes[1].points.begin(), reordered.holes[1].points.end());

    const CanonicalGeometry first_geometry = canonical_geometry(first);
    const CanonicalGeometry reordered_geometry = canonical_geometry(reordered);
    REQUIRE(first_geometry == reordered_geometry);
    REQUIRE(canonical_hash(first_geometry) == canonical_hash(reordered_geometry));

    reordered.holes[0].points[0].x() += 1;
    const CanonicalGeometry changed_geometry = canonical_geometry(reordered);
    REQUIRE_FALSE(first_geometry == changed_geometry);

    // An intentionally shared hash bucket must still select by exact geometry.
    LayerCache collision_bucket;
    collision_bucket[42].push_back({ first_geometry, {}, 7 });
    collision_bucket[42].push_back({ changed_geometry, {}, 9 });
    REQUIRE(find_cached_layer(collision_bucket, 42, first_geometry)->source_layer == 7);
    REQUIRE(find_cached_layer(collision_bucket, 42, changed_geometry)->source_layer == 9);
    CanonicalGeometry absent_geometry = first_geometry;
    absent_geometry.outer.front().first += 1;
    REQUIRE(find_cached_layer(collision_bucket, 42, absent_geometry) == nullptr);
}

TEST_CASE("Continuous Fermat internet model corpus", "[.][ContinuousFermat][internet-corpus]")
{
    const std::string manifest_value = env_string("CONTINUOUS_FERMAT_CORPUS_MANIFEST");
    if (manifest_value.empty())
        SKIP("set CONTINUOUS_FERMAT_CORPUS_MANIFEST to run the external model corpus");

    const stdfs::path manifest = stdfs::absolute(stdfs::u8path(manifest_value));
    const double layer_height = env_double("CONTINUOUS_FERMAT_CORPUS_LAYER_HEIGHT", 0.20);
    const double line_width = env_double("CONTINUOUS_FERMAT_CORPUS_LINE_WIDTH", 0.40);
    const double nozzle_diameter = env_double("CONTINUOUS_FERMAT_CORPUS_NOZZLE", 0.40);
    const double max_line_width = env_double("CONTINUOUS_FERMAT_CORPUS_MAX_LINE_WIDTH", 0.80);
    const double max_dimension = env_double("CONTINUOUS_FERMAT_CORPUS_MAX_DIMENSION_MM", 0.0);
    const size_t max_layers = env_size("CONTINUOUS_FERMAT_CORPUS_MAX_LAYERS", 2000);
    const size_t max_models = env_size("CONTINUOUS_FERMAT_CORPUS_MAX_MODELS", std::numeric_limits<size_t>::max());
    const size_t shard_index = env_size("CONTINUOUS_FERMAT_CORPUS_SHARD_INDEX", 0);
    const size_t shard_count = std::max<size_t>(1, env_size("CONTINUOUS_FERMAT_CORPUS_SHARD_COUNT", 1));
    const size_t progress_every = env_size("CONTINUOUS_FERMAT_CORPUS_PROGRESS_EVERY", 1);
    const double layer_budget_ms = env_double("CONTINUOUS_FERMAT_CORPUS_LAYER_BUDGET_MS", 0.0);
    const double model_budget_ms = env_double("CONTINUOUS_FERMAT_CORPUS_MODEL_BUDGET_MS", 0.0);
    const bool deduplicate_layers = env_bool("CONTINUOUS_FERMAT_CORPUS_DEDUPLICATE_LAYERS", true);
    const bool require_all = env_bool("CONTINUOUS_FERMAT_CORPUS_REQUIRE_ALL", false);
    REQUIRE(layer_height > 0.0);
    REQUIRE(line_width > 0.0);
    REQUIRE(nozzle_diameter > 0.0);
    REQUIRE(max_line_width >= line_width);
    REQUIRE(shard_index < shard_count);
    REQUIRE(layer_budget_ms >= 0.0);
    REQUIRE(model_budget_ms >= 0.0);

    std::vector<CorpusEntry> entries = read_manifest(manifest);
    const stdfs::path output_prefix = stdfs::u8path(env_string(
        "CONTINUOUS_FERMAT_CORPUS_OUTPUT", (manifest.parent_path() / (manifest.stem().u8string() + "-results")).u8string()));
    if (!output_prefix.parent_path().empty())
        stdfs::create_directories(output_prefix.parent_path());
    std::ofstream jsonl(output_prefix.string() + ".jsonl", std::ios::trunc);
    std::ofstream csv(output_prefix.string() + ".csv", std::ios::trunc);
    REQUIRE(jsonl.good());
    REQUIRE(csv.good());
    csv << "record,id,path,layer,z_mm,islands,holes,area_mm2,overlap_mm2,eligible,attempted,passed,quality_ok,category,reason,generation_ms,validation_ms,path_points,closed,exact_coverage,coverage,outside,material,redeposition,containment,crossings,spacing,bead_overlaps,turnbacks,thin_feature_class,thin_feature_threshold_mm,load_ms,slice_ms,total_ms,effective_scale,source_url,layer_elapsed_ms,cache_hit,reused_from_layer,unique_geometries,cache_hits\n";

    size_t processed_models = 0;
    size_t eligible_models = 0;
    size_t passed_models = 0;
    size_t failed_layers = 0;
    size_t ineligible_models = 0;
    size_t harness_errors = 0;
    size_t total_unique_geometries = 0;
    size_t total_cache_hits = 0;
    for (size_t manifest_index = 0; manifest_index < entries.size() && processed_models < max_models; ++manifest_index) {
        if (manifest_index % shard_count != shard_index)
            continue;
        CorpusEntry &entry = entries[manifest_index];
        ++processed_models;
        const Clock::time_point model_start = Clock::now();
        double load_ms = 0.0;
        double slice_ms = 0.0;
        double effective_scale = entry.scale;
        std::string model_category;
        std::string model_reason;
        std::vector<LayerResult> results;
        std::vector<bool> layer_written;
        bool model_eligible = false;
        bool model_passed = false;
        bool model_quality_ok = false;
        size_t unique_geometries = 0;
        size_t cache_hits = 0;

        jsonl << "{\"record\":\"model_start\",\"id\":\"" << json_escape(entry.id) << "\",\"path\":\""
              << json_escape(entry.path.u8string()) << "\",\"model_number\":" << processed_models << "}\n";
        jsonl.flush();
        std::cout << "corpus model " << processed_models << ": " << entry.id << " started" << std::endl;

        const auto emit_layer = [&](const size_t index) {
            if (index >= results.size() || (index < layer_written.size() && layer_written[index]))
                return;
            write_layer_json(jsonl, entry, results[index]);
            write_layer_csv(csv, entry, results[index]);
            jsonl.flush();
            csv.flush();
            if (index < layer_written.size())
                layer_written[index] = true;
            if (progress_every != 0 && ((index + 1) % progress_every == 0 || index + 1 == results.size())) {
                const LayerResult &result = results[index];
                const char *work = result.cache_hit ? "cached" : result.attempted ? "unique" : "skipped";
                std::cout << "corpus model " << processed_models << " " << entry.id << " layer " << (index + 1) << '/' << results.size()
                          << ' ' << work << " category=" << result.category
                          << " layer_ms=" << std::fixed << std::setprecision(1) << result.layer_elapsed_ms
                          << " model_ms=" << elapsed_ms(model_start) << std::defaultfloat << std::endl;
            }
        };

        try {
            const Clock::time_point load_start = Clock::now();
            Model model = Model::read_from_file(entry.path.u8string(), nullptr, nullptr,
                LoadStrategy::AddDefaultInstances | LoadStrategy::Silence);
            TriangleMesh mesh = model.mesh();
            if (mesh.empty())
                throw std::runtime_error("loaded mesh is empty");
            if (entry.scale <= 0.0 || !std::isfinite(entry.scale))
                throw std::runtime_error("manifest scale must be finite and positive");
            if (entry.scale != 1.0)
                mesh.scale(float(entry.scale));
            BoundingBoxf3 bbox = mesh.bounding_box();
            if (!bbox.defined || !bbox.min.allFinite() || !bbox.max.allFinite())
                throw std::runtime_error("mesh bounding box is invalid");
            if (max_dimension > 0.0) {
                const double largest = bbox.size().maxCoeff();
                if (largest > max_dimension) {
                    const double fit_scale = max_dimension / largest;
                    mesh.scale(float(fit_scale));
                    effective_scale *= fit_scale;
                    bbox = mesh.bounding_box();
                }
            }
            mesh.translate(0.0f, 0.0f, float(-bbox.min.z()));
            bbox = mesh.bounding_box();
            load_ms = elapsed_ms(load_start);

            const double height = bbox.size().z();
            if (!(height > 0.0))
                throw std::runtime_error("mesh has zero height");
            const size_t layer_count = size_t(std::ceil(height / layer_height - 1e-9));
            if (layer_count == 0 || layer_count > max_layers) {
                model_category = "layer_limit";
                model_reason = "layer count " + std::to_string(layer_count) + " exceeds limit " + std::to_string(max_layers);
                ++ineligible_models;
                ++harness_errors;
            } else {
                std::vector<double> zs;
                zs.reserve(layer_count);
                for (size_t i = 0; i < layer_count; ++i) {
                    const double bottom = double(i) * layer_height;
                    const double top = std::min(height, bottom + layer_height);
                    zs.emplace_back(0.5 * (bottom + top));
                }
                const Clock::time_point slice_start = Clock::now();
                const std::vector<ExPolygons> slices = mesh.slice(zs);
                slice_ms = elapsed_ms(slice_start);
                results.resize(slices.size());
                layer_written.assign(slices.size(), false);
                model_eligible = !slices.empty();
                ExPolygons previous;
                for (size_t i = 0; i < slices.size(); ++i) {
                    LayerResult &result = results[i];
                    result.index = i;
                    result.z = zs[i];
                    result.islands = slices[i].size();
                    result.area_mm2 = area_mm2(slices[i]);
                    for (const ExPolygon &island : slices[i])
                        result.holes += island.holes.size();
                    if (slices[i].size() != 1 || result.area_mm2 <= 0.0) {
                        result.category = slices[i].empty() ? "empty_layer" : "multiple_islands";
                        result.reason = slices[i].empty() ? "nominal layer has no printable area" : "layer does not contain exactly one island";
                        model_eligible = false;
                    } else if (!previous.empty()) {
                        const ExPolygons overlap = intersection_ex(previous, slices[i]);
                        result.overlap_mm2 = area_mm2(overlap);
                        if (result.overlap_mm2 <= 0.0) {
                            result.category = "no_positive_overlap";
                            result.reason = "layer has no positive-area overlap with the previous layer";
                            model_eligible = false;
                        }
                    }
                    result.eligible = result.category.empty();
                    previous = slices[i];
                }

                if (model_eligible) {
                    ++eligible_models;
                    model_passed = true;
                    model_quality_ok = true;
                    const Flow flow { float(line_width), float(layer_height), float(nozzle_diameter) };
                    LayerCache cache;
                    bool stop_for_model_budget = false;
                    bool budget_exceeded = false;
                    for (size_t i = 0; i < slices.size(); ++i) {
                        LayerResult &result = results[i];
                        const Clock::time_point layer_start = Clock::now();
                        if (stop_for_model_budget || (model_budget_ms > 0.0 && elapsed_ms(model_start) >= model_budget_ms)) {
                            stop_for_model_budget = true;
                            budget_exceeded = true;
                            result.category = "model_budget_exceeded";
                            result.reason = "model elapsed budget was exhausted before this layer";
                            result.passed = false;
                        } else {
                            CanonicalGeometry geometry;
                            uint64_t geometry_hash = 0;
                            const CachedLayerResult *cached = nullptr;
                            if (deduplicate_layers) {
                                geometry = canonical_geometry(slices[i].front());
                                geometry_hash = canonical_hash(geometry);
                                cached = find_cached_layer(cache, geometry_hash, geometry);
                            }

                            if (cached != nullptr) {
                                const size_t index = result.index;
                                const double z = result.z;
                                const size_t islands = result.islands;
                                const size_t holes = result.holes;
                                const double layer_area = result.area_mm2;
                                const double overlap_area = result.overlap_mm2;
                                const bool eligible = result.eligible;
                                result = cached->result;
                                result.index = index;
                                result.z = z;
                                result.islands = islands;
                                result.holes = holes;
                                result.area_mm2 = layer_area;
                                result.overlap_mm2 = overlap_area;
                                result.eligible = eligible;
                                result.attempted = false;
                                result.cache_hit = true;
                                result.reused_from_layer = cached->source_layer;
                                result.generation_ms = 0.0;
                                result.validation_ms = 0.0;
                                ++cache_hits;
                            } else {
                                result.attempted = true;
                                Clock::time_point operation_start = Clock::now();
                                bool validating = false;
                                try {
                                    const ContinuousFermat::GeneratedPath generated =
                                        ContinuousFermat::generate_layer_path_with_metadata(slices[i], flow, max_line_width);
                                    result.generation_ms = elapsed_ms(operation_start);
                                    result.path_points = generated.path.points.size();
                                    validating = true;
                                    operation_start = Clock::now();
                                    result.validation = ContinuousFermat::validate_layer_path(
                                        slices[i], flow, generated.path, generated.extrusion_multipliers, max_line_width);
                                    result.validation_ms = elapsed_ms(operation_start);
                                    result.passed = result.validation.emittable;
                                    result.quality_ok = result.validation.ok;
                                    result.reason = result.validation.reason;
                                    result.category = result.quality_ok ? "pass" :
                                        result.passed ? "advisory_" + validation_category(result.validation) :
                                                        validation_category(result.validation);
                                } catch (const std::exception &e) {
                                    if (validating)
                                        result.validation_ms = elapsed_ms(operation_start);
                                    else
                                        result.generation_ms = elapsed_ms(operation_start);
                                    result.category = "layer_exception";
                                    result.reason = e.what();
                                    result.passed = false;
                                } catch (...) {
                                    if (validating)
                                        result.validation_ms = elapsed_ms(operation_start);
                                    else
                                        result.generation_ms = elapsed_ms(operation_start);
                                    result.category = "layer_exception";
                                    result.reason = "unknown exception";
                                    result.passed = false;
                                }
                                result.layer_elapsed_ms = elapsed_ms(layer_start);
                                if (layer_budget_ms > 0.0 && result.layer_elapsed_ms > layer_budget_ms) {
                                    budget_exceeded = true;
                                    result.category = "layer_budget_exceeded";
                                    result.reason = "layer elapsed time exceeded the configured discovery budget";
                                    result.passed = false;
                                }
                                ++unique_geometries;
                                if (deduplicate_layers)
                                    cache[geometry_hash].push_back({ std::move(geometry), result, i });
                            }

                            result.layer_elapsed_ms = elapsed_ms(layer_start);
                            if (layer_budget_ms > 0.0 && result.layer_elapsed_ms > layer_budget_ms && result.category != "layer_budget_exceeded") {
                                budget_exceeded = true;
                                result.category = "layer_budget_exceeded";
                                result.reason = "layer elapsed time exceeded the configured discovery budget";
                                result.passed = false;
                            }
                            if (model_budget_ms > 0.0 && elapsed_ms(model_start) > model_budget_ms) {
                                stop_for_model_budget = true;
                                budget_exceeded = true;
                                result.category = "model_budget_exceeded";
                                result.reason = "model elapsed time exceeded the configured discovery budget";
                                result.passed = false;
                            }
                        }
                        result.layer_elapsed_ms = elapsed_ms(layer_start);
                        if (!result.passed) {
                            ++failed_layers;
                            model_passed = false;
                        }
                        if (!result.quality_ok)
                            model_quality_ok = false;
                        emit_layer(i);
                    }
                    if (model_passed) {
                        ++passed_models;
                        model_category = model_quality_ok ? "pass" : "pass_with_advisories";
                    }
                    else if (budget_exceeded) {
                        model_category = "budget_exceeded";
                        model_reason = "one or more layers exceeded a discovery elapsed-time budget";
                    }
                    else {
                        model_category = "fermat_failure";
                        model_reason = "one or more eligible layers failed generation or validation";
                    }
                } else {
                    ++ineligible_models;
                    model_category = "topology_ineligible";
                    model_reason = "one-island/positive-overlap eligibility rule failed";
                }
            }
        } catch (const std::exception &e) {
            model_category = "load_or_slice_exception";
            model_reason = e.what();
            ++ineligible_models;
            ++harness_errors;
        } catch (...) {
            model_category = "load_or_slice_exception";
            model_reason = "unknown exception";
            ++ineligible_models;
            ++harness_errors;
        }

        for (size_t i = 0; i < results.size(); ++i)
            emit_layer(i);
        total_unique_geometries += unique_geometries;
        total_cache_hits += cache_hits;
        const double model_ms = elapsed_ms(model_start);
        jsonl << std::setprecision(10) << "{\"record\":\"model\",\"id\":\"" << json_escape(entry.id)
              << "\",\"path\":\"" << json_escape(entry.path.u8string()) << "\",\"source_url\":\"" << json_escape(entry.source_url)
              << "\",\"scale\":" << effective_scale << ",\"layers\":" << results.size()
              << ",\"eligible\":" << (model_eligible ? "true" : "false") << ",\"passed\":" << (model_passed ? "true" : "false")
              << ",\"quality_ok\":" << (model_quality_ok ? "true" : "false")
              << ",\"category\":\"" << json_escape(model_category) << "\",\"reason\":\"" << json_escape(model_reason)
              << "\",\"load_ms\":" << load_ms << ",\"slice_ms\":" << slice_ms << ",\"total_ms\":" << model_ms
              << ",\"unique_geometries\":" << unique_geometries << ",\"cache_hits\":" << cache_hits << "}\n";
        csv << std::setprecision(10) << "model," << csv_escape(entry.id) << ',' << csv_escape(entry.path.u8string());
        for (size_t i = 0; i < 7; ++i)
            csv << ',';
        csv << model_eligible << ",," << model_passed << ',' << model_quality_ok << ','
            << csv_escape(model_category) << ',' << csv_escape(model_reason);
        for (size_t i = 0; i < 17; ++i)
            csv << ',';
        csv << load_ms << ',' << slice_ms << ',' << model_ms << ',' << effective_scale << ',' << csv_escape(entry.source_url) << ",,,,"
            << unique_geometries << ',' << cache_hits << '\n';
        jsonl.flush();
        csv.flush();
    }

    jsonl << "{\"record\":\"run\",\"processed_models\":" << processed_models << ",\"eligible_models\":" << eligible_models
          << ",\"passed_models\":" << passed_models << ",\"ineligible_models\":" << ineligible_models
          << ",\"harness_errors\":" << harness_errors << ",\"failed_layers\":" << failed_layers << ",\"shard_index\":" << shard_index
          << ",\"shard_count\":" << shard_count << ",\"unique_geometries\":" << total_unique_geometries
          << ",\"cache_hits\":" << total_cache_hits << "}\n";
    INFO("results: " << output_prefix.u8string() << ".jsonl and .csv");
    INFO("processed=" << processed_models << " eligible=" << eligible_models << " passed=" << passed_models
         << " ineligible=" << ineligible_models << " harness_errors=" << harness_errors << " failed_layers=" << failed_layers);
    REQUIRE(processed_models > 0);
    if (require_all) {
        CHECK(harness_errors == 0);
        CHECK(failed_layers == 0);
        CHECK(passed_models == eligible_models);
    }
}
