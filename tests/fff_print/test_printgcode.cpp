#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
	#include <Windows.h>
#endif

#include <catch2/catch_all.hpp>

#include "libslic3r/libslic3r.h"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Layer.hpp"

#include "test_data.hpp"

#include <algorithm>
#include <boost/regex.hpp>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

using namespace Slic3r;
using namespace Slic3r::Test;

boost::regex perimeters_regex("G1 X[-0-9.]* Y[-0-9.]* E[-0-9.]* ; perimeter");
boost::regex infill_regex("G1 X[-0-9.]* Y[-0-9.]* E[-0-9.]* ; infill");
boost::regex skirt_regex("G1 X[-0-9.]* Y[-0-9.]* E[-0-9.]* ; skirt");

namespace {

std::map<char, double> parse_numeric_words(const std::string &code)
{
    std::istringstream input(code);
    std::string token;
    input >> token; // command
    std::map<char, double> words;
    while (input >> token) {
        if (token.size() < 2 || !std::isalpha(static_cast<unsigned char>(token.front())))
            continue;
        size_t consumed = 0;
        const double value = std::stod(token.substr(1), &consumed);
        if (consumed != token.size() - 1)
            throw std::runtime_error("Malformed numeric G-code word: " + token);
        words.emplace(char(std::toupper(static_cast<unsigned char>(token.front()))), value);
    }
    return words;
}

DynamicPrintConfig strict_continuous_export_config()
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"spiral_mode", false},
        {"spiral_hybrid_non_crossing", true},
        {"gcode_flavor", "marlin2"},
        {"print_sequence", "by layer"},
        {"nozzle_diameter", "1.2"},
        {"outer_wall_line_width", 1.2},
        {"initial_layer_line_width", 1.2},
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"skirt_loops", 0},
        {"skirt_height", 0},
        {"brim_type", "no_brim"},
        {"brim_width", 0.0},
        {"enable_prime_tower", false},
        {"exclude_object", false},
        {"enable_power_loss_recovery", "disable"},
        {"enable_support", false},
        {"raft_layers", 0},
        {"scan_first_layer", false},
        {"ironing_type", "no ironing"},
        {"before_layer_change_gcode", ""},
        {"layer_change_gcode", ""},
        {"time_lapse_gcode", ""},
        {"change_extrusion_role_gcode", ""},
        {"file_start_gcode", ""},
        {"machine_start_gcode", ""},
        {"machine_end_gcode", ""},
        {"filament_start_gcode", ""},
        {"filament_end_gcode", ""},
        {"nozzle_temperature_initial_layer", "200"},
        {"nozzle_temperature", "200"},
        {"use_relative_e_distances", true},
        {"retract_when_changing_layer", "1"},
        {"retraction_length", "2"},
        {"z_hop", "1"},
        {"wipe", "1"},
        {"enable_arc_fitting", true},
        {"resonance_avoidance", true},
        {"print_flow_ratio", 1.30},
        {"outer_wall_flow_ratio", 1.40},
        {"first_layer_flow_ratio", 1.40},
        {"filament_flow_ratio", "1"},
        {"filament_max_volumetric_speed", "2"},
        {"filament_adaptive_volumetric_speed", "0"},
        {"max_volumetric_extrusion_rate_slope", 0.0},
        {"enable_pressure_advance", "0"},
        {"adaptive_pressure_advance", "0"},
        {"slow_down_for_layer_cooling", "1"},
        {"slow_down_layer_time", "1000"},
        {"emit_machine_limits_to_gcode", false},
        {"gcode_add_line_number", false},
    });
    return config;
}

void apply_continuous_prism(
    Print &print,
    Model &model,
    const DynamicPrintConfig &config,
    const Vec3d &offset = Vec3d(20.0, 20.0, 0.0),
    bool allow_development_export_for_test = true,
    float xy_scale = 1.2f,
    float z_scale = 0.03f)
{
    TriangleMesh prism = Test::mesh(TestMesh::cube_20x20x20);
    prism.scale(Vec3f(xy_scale, xy_scale, z_scale));
    ModelObject *object = model.add_object();
    object->name = "continuous-test-prism.stl";
    object->add_volume(std::move(prism));
    ModelInstance *instance = object->add_instance();
    instance->set_offset(offset);
    object->ensure_on_bed();
    print.auto_assign_extruders(object);
    print.apply(model, config);
    if (allow_development_export_for_test)
        print.enable_continuous_slicing_development_export_for_tests();
}

bool continuous_export_throws_slicing_error(Print &print)
{
    const std::string path = boost::filesystem::unique_path().string();
    const auto cleanup = [&path]() {
        std::remove(path.c_str());
        std::remove((path + ".tmp").c_str());
    };
    try {
        print.set_status_silent();
        print.process();
        print.export_gcode(path, nullptr, nullptr);
    } catch (const Slic3r::SlicingError &) {
        cleanup();
        return true;
    } catch (...) {
        cleanup();
        throw;
    }
    cleanup();
    return false;
}

} // namespace

TEST_CASE("Continuous slicing permits validated printer-ready G-code export", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    Print print;
    Model model;
    apply_continuous_prism(print, model, config, Vec3d(20.0, 20.0, 0.0), false);
    REQUIRE(print.validate().string.empty());
    print.set_status_silent();
    print.process();

    const std::string path = boost::filesystem::unique_path().string();
    REQUIRE_NOTHROW(print.export_gcode(path, nullptr, nullptr));
    {
        std::ifstream input(path, std::ios::binary);
        const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
        CHECK(contents.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
        CHECK(contents.find(";_CONTINUOUS_FERMAT_END") != std::string::npos);
    }
    CHECK_FALSE(boost::filesystem::exists(path + ".tmp"));
    std::remove(path.c_str());

    const std::string preview_path = boost::filesystem::unique_path().string();
    REQUIRE_NOTHROW(print.export_gcode(
        preview_path,
        nullptr,
        nullptr,
        GCodeExportPurpose::ContinuousResearchPreview));
    {
        std::ifstream input(preview_path, std::ios::binary);
        const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
        CHECK(contents.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
        CHECK(contents.find("; spiral_mode = 0") != std::string::npos);
        CHECK(contents.find("; spiral_hybrid_non_crossing = 1") != std::string::npos);
    }
    std::remove(preview_path.c_str());
}

TEST_CASE("Continuous slicing exports structurally valid paths with geometric warnings", "[PrintGCode][continuous][advisory]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    Print print;
    Model model;
    // The 6 mm square is a known closed, finite candidate whose physical
    // footprint misses at least one fixed geometric quality target at a
    // 1.2 mm line width.
    apply_continuous_prism(print, model, config, Vec3d(20.0, 20.0, 0.0), false, 0.3f, 0.03f);
    REQUIRE(print.validate().string.empty());
    print.set_status_silent();
    REQUIRE_NOTHROW(print.process());

    const std::string path = boost::filesystem::unique_path().string();
    REQUIRE_NOTHROW(print.export_gcode(path, nullptr, nullptr));
    std::ifstream input(path, std::ios::binary);
    const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    CHECK(contents.find(";_CONTINUOUS_FERMAT_VALIDATION_WARNING") != std::string::npos);
    CHECK(contents.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
    CHECK(contents.find(";_CONTINUOUS_FERMAT_END") != std::string::npos);
    std::remove(path.c_str());
}

TEST_CASE("Continuous slicing ignores a wider first-layer line on a 0.4 mm nozzle", "[PrintGCode][continuous][regression]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set_deserialize_strict("nozzle_diameter", "0.4");
    config.set("outer_wall_line_width", 0.4);
    config.set("initial_layer_line_width", 0.5);

    Print print;
    Model model;
    apply_continuous_prism(print, model, config, Vec3d(20.0, 20.0, 0.0), true, 0.6f, 0.5f);
    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_END") != std::string::npos);
}

TEST_CASE("Continuous slicing processes a hole-free cylindrical prism", "[PrintGCode][continuous][regression]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    Print print;
    Model model;
    ModelObject *object = model.add_object();
    object->name = "continuous-test-cylinder.stl";
    object->add_volume(make_cylinder(24.0, 0.6, PI / 48.0));
    ModelInstance *instance = object->add_instance();
    instance->set_offset(Vec3d(50.0, 50.0, 0.0));
    object->ensure_on_bed();
    print.auto_assign_extruders(object);
    print.apply(model, config);
    print.enable_continuous_slicing_development_export_for_tests();

    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
    REQUIRE(std::count(gcode.begin(), gcode.end(), '\n') > 3);
}

TEST_CASE("Continuous slicing processes an overlapping pyramid cross-section", "[PrintGCode][continuous][regression]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    Print print;
    Model model;
    ModelObject *object = model.add_object();
    object->name = "continuous-test-pyramid.stl";
    object->add_volume(Test::mesh(TestMesh::pyramid));
    ModelInstance *instance = object->add_instance();
    instance->set_offset(Vec3d(50.0, 50.0, 0.0));
    object->ensure_on_bed();
    print.auto_assign_extruders(object);
    print.apply(model, config);
    print.enable_continuous_slicing_development_export_for_tests();

    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_BEGIN") != std::string::npos);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_END") != std::string::npos);
}

TEST_CASE("Cached Continuous Fermat artifacts remain blocked after the mode is disabled", "[PrintGCode][continuous]")
{
    Print print;

    const std::string path = boost::filesystem::unique_path().string();
    {
        std::ofstream output(path, std::ios::binary);
        REQUIRE(output.is_open());
        output << ";_CONTINUOUS_FERMAT_BEGIN\nG1 X1 Y1 E1\n;_CONTINUOUS_FERMAT_END\n";
    }
    CHECK_NOTHROW(print.throw_if_continuous_slicing_artifact_blocked(path, true));
    REQUIRE_THROWS_WITH(
        print.throw_if_continuous_slicing_artifact_blocked(path),
        Catch::Matchers::ContainsSubstring("cached/imported Continuous Fermat G-code is intentionally disabled"));
    std::remove(path.c_str());
}

TEST_CASE("Continuous slicing export preserves the serialized section contract", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    // The G-code footer serializes the complete configuration.  Generic
    // vector-enum defaults must retain their key map when copied through the
    // static configuration cache.
    REQUIRE(config.opt_serialize("extruder_type") == "Direct Drive");

    Print print;
    Model model;
    print.set_plate_origin(Vec3d(0.0006, 0.0004, 0.0));
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);
    if (const char *artifact_path = std::getenv("ORCA_CONTINUOUS_TEST_GCODE");
        artifact_path != nullptr && *artifact_path != '\0') {
        std::ofstream artifact(artifact_path, std::ios::binary);
        REQUIRE(artifact.good());
        artifact.write(gcode.data(), std::streamsize(gcode.size()));
        REQUIRE(artifact.good());
    }

    size_t begin_count = 0;
    size_t end_count = 0;
    size_t layer_count = 0;
    bool inside = false;
    bool after_section = false;
    bool saw_z_transition = false;
    bool saw_forbidden_transition_command = false;
    double current_x = std::numeric_limits<double>::quiet_NaN();
    double current_y = std::numeric_limits<double>::quiet_NaN();
    double current_z = std::numeric_limits<double>::quiet_NaN();
    double current_f = std::numeric_limits<double>::quiet_NaN();
    double section_start_x = 0.0;
    double section_start_y = 0.0;
    double section_start_z = 0.0;
    double section_volume = 0.0;
    size_t section_extrusion_moves = 0;
    const double filament_area = PI * 1.75 * 1.75 * 0.25;
    const double expected_layer_volume = 24.0 * 24.0 * 0.2;
    const double max_volumetric_speed = 2.0;
    std::istringstream stream(gcode);
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line == ";LAYER_CHANGE" || line == "; CHANGE_LAYER")
            ++layer_count;
        if (line == ";_CONTINUOUS_FERMAT_BEGIN") {
            if (after_section) {
                CHECK(saw_z_transition);
                CHECK_FALSE(saw_forbidden_transition_command);
            }
            REQUIRE(std::isfinite(current_x));
            REQUIRE(std::isfinite(current_y));
            REQUIRE(std::isfinite(current_z));
            REQUIRE(current_z > 0.0);
            section_start_x = current_x;
            section_start_y = current_y;
            section_start_z = current_z;
            section_volume = 0.0;
            section_extrusion_moves = 0;
            inside = true;
            after_section = false;
            saw_z_transition = false;
            saw_forbidden_transition_command = false;
            ++begin_count;
            continue;
        }
        if (line == ";_CONTINUOUS_FERMAT_END") {
            REQUIRE(inside);
            REQUIRE(section_extrusion_moves > 0);
            CHECK(current_x == section_start_x);
            CHECK(current_y == section_start_y);
            CHECK(current_z == section_start_z);
            CHECK(section_volume / expected_layer_volume >= 0.980);
            CHECK(section_volume / expected_layer_volume <= 1.020);
            inside = false;
            after_section = true;
            ++end_count;
            continue;
        }

        const std::string code = line.substr(0, line.find(';'));
        const bool g0 = code.rfind("G0 ", 0) == 0;
        const bool g1 = code.rfind("G1 ", 0) == 0;
        const bool has_x = code.find(" X") != std::string::npos;
        const bool has_y = code.find(" Y") != std::string::npos;
        const bool has_z = code.find(" Z") != std::string::npos;
        const bool has_e = code.find(" E") != std::string::npos;
        const std::map<char, double> words = (g0 || g1) ? parse_numeric_words(code) : std::map<char, double>{};
        if (auto it = words.find('F'); it != words.end()) {
            REQUIRE(std::isfinite(it->second));
            REQUIRE(it->second > 0.0);
            current_f = it->second;
        }
        const double next_x = words.count('X') ? words.at('X') : current_x;
        const double next_y = words.count('Y') ? words.at('Y') : current_y;
        const double next_z = words.count('Z') ? words.at('Z') : current_z;
        if (inside) {
            CHECK_FALSE(g0);
            CHECK(code.rfind("G2 ", 0) != 0);
            CHECK(code.rfind("G3 ", 0) != 0);
            CHECK_FALSE(has_z);
            if (g1 && (has_x || has_y)) {
                REQUIRE(has_e);
                REQUIRE(words.count('E') == 1);
                REQUIRE(std::isfinite(words.at('E')));
                REQUIRE(words.at('E') > 0.0);
                REQUIRE(std::isfinite(current_x));
                REQUIRE(std::isfinite(current_y));
                REQUIRE(std::isfinite(next_x));
                REQUIRE(std::isfinite(next_y));
                REQUIRE(std::isfinite(current_f));
                const double length = std::hypot(next_x - current_x, next_y - current_y);
                INFO("section=" << begin_count << " line='" << line << "' start=" << current_x << "," << current_y <<
                     " end=" << next_x << "," << next_y << " F=" << current_f);
                REQUIRE(length > 0.0);
                const double volume = words.at('E') * filament_area;
                const double volumetric_speed = volume * current_f / (60.0 * length);
                CHECK(volumetric_speed <= max_volumetric_speed + 1e-9);
                section_volume += volume;
                ++section_extrusion_moves;
            } else if (g1) {
                CHECK_FALSE(has_e);
            }
        } else if (after_section) {
            if ((g0 || g1) && has_z && !has_x && !has_y && !has_e)
                saw_z_transition = true;
            const bool forbidden_motion = (g0 || g1) && (has_x || has_y || has_e);
            saw_forbidden_transition_command |= forbidden_motion ||
                code.rfind("G10", 0) == 0 || code.rfind("G11", 0) == 0 || code.rfind("G92", 0) == 0;
        }
        if (g0 || g1) {
            current_x = next_x;
            current_y = next_y;
            current_z = next_z;
        }
    }

    REQUIRE_FALSE(inside);
    REQUIRE(begin_count == end_count);
    REQUIRE(begin_count == layer_count);
    REQUIRE(begin_count == 3);
    if (after_section)
        CHECK_FALSE(saw_forbidden_transition_command);
    REQUIRE(gcode.find("G21 ; continuous slicing: millimetres") != std::string::npos);
    REQUIRE(gcode.find("G90 ; continuous slicing: absolute XYZ") != std::string::npos);
    REQUIRE(gcode.find("M83 ; continuous slicing: relative E") != std::string::npos);
    REQUIRE(gcode.find("T0 ; continuous slicing: select the certified tool") != std::string::npos);
    REQUIRE(gcode.find("M200 D0 ; continuous slicing: use filament-length E") != std::string::npos);
    REQUIRE(gcode.find("M220 S100 ; continuous slicing: reset speed override") != std::string::npos);
    REQUIRE(gcode.find("M221 S100 ; continuous slicing: reset flow override") != std::string::npos);
    REQUIRE(gcode.find("M900 K0 ; continuous slicing: disable linear advance") != std::string::npos);
    REQUIRE(gcode.find("M413 S0 ; continuous slicing: disable power-loss recovery") != std::string::npos);
    REQUIRE(gcode.find("G28 ; continuous slicing: home all axes") != std::string::npos);
    REQUIRE(gcode.find("M109 R200 ; continuous slicing: wait for certified nozzle temperature while heating or cooling") != std::string::npos);

    const size_t m900_position = gcode.find("M900 K0 ; continuous slicing: disable linear advance");
    const size_t home_position = gcode.find("G28 ; continuous slicing: home all axes");
    const size_t wait_position = gcode.find("M109 R200 ; continuous slicing: wait for certified nozzle temperature while heating or cooling");
    const size_t begin_position = gcode.find(";_CONTINUOUS_FERMAT_BEGIN");
    REQUIRE(m900_position < home_position);
    REQUIRE(home_position < wait_position);
    REQUIRE(wait_position < begin_position);
}

TEST_CASE("Continuous slicing emits a Klipper-native protected preview contract", "[PrintGCode][continuous][klipper]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set_deserialize_strict("gcode_flavor", "klipper");
    config.set("file_start_gcode", "PRINT_START");
    config.set("machine_start_gcode", "START_PRINT BED=60 HOTEND=200");
    config.set("machine_end_gcode", "END_PRINT");

    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);
    if (const char *artifact_path = std::getenv("ORCA_CONTINUOUS_KLIPPER_TEST_GCODE");
        artifact_path != nullptr && *artifact_path != '\0') {
        std::ofstream artifact(artifact_path, std::ios::binary);
        REQUIRE(artifact.good());
        artifact.write(gcode.data(), std::streamsize(gcode.size()));
        REQUIRE(artifact.good());
    }

    const size_t home = gcode.find("G28 ; continuous slicing: home all axes");
    const size_t hotend_wait = gcode.find("M109 S200 ; continuous slicing: wait for certified nozzle temperature");
    const size_t first_section = gcode.find(";_CONTINUOUS_FERMAT_BEGIN");
    REQUIRE(home != std::string::npos);
    REQUIRE(hotend_wait != std::string::npos);
    REQUIRE(first_section != std::string::npos);
    CHECK(home < hotend_wait);
    CHECK(hotend_wait < first_section);
    CHECK(gcode.find("SET_PRESSURE_ADVANCE ADVANCE=0 ; continuous slicing: disable pressure advance") != std::string::npos);
    CHECK(gcode.find("M220 S100 ; continuous slicing: reset speed override") != std::string::npos);
    CHECK(gcode.find("M221 S100 ; continuous slicing: reset flow override") != std::string::npos);
    CHECK((gcode.find("M190 S") != std::string::npos || gcode.find("M140 S0 ; continuous slicing: heated bed disabled") != std::string::npos));
    CHECK(gcode.find("\nPRINT_START\n") == std::string::npos);
    CHECK(gcode.find("\nSTART_PRINT BED=60 HOTEND=200\n") == std::string::npos);
    CHECK(gcode.find("\nEND_PRINT\n") == std::string::npos);
    CHECK(gcode.find("\nG21") == std::string::npos);
    CHECK(gcode.find("\nT0 ") == std::string::npos);
    CHECK(gcode.find("\nM200 ") == std::string::npos);
    CHECK(gcode.find("\nM900 ") == std::string::npos);
    CHECK(gcode.find("\nM413 ") == std::string::npos);
    CHECK(gcode.find("ACCEL_TO_DECEL") == std::string::npos);
    CHECK(gcode.find("\nM73 ") == std::string::npos);
}

TEST_CASE("Continuous slicing preserves filament flow calibration and omits incompatible features", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set_deserialize_strict("gcode_flavor", "klipper");
    config.set_deserialize_strict("filament_flow_ratio", "0.93");
    config.set_deserialize_strict("brim_type", "outer_only");
    config.set("brim_width", 3.0);
    config.set("exclude_object", true);
    config.set_deserialize_strict("enable_power_loss_recovery", "enable");
    config.set("emit_machine_limits_to_gcode", true);
    config.set("gcode_add_line_number", true);
    config.set("max_volumetric_extrusion_rate_slope", 1.0);
    config.set_deserialize_strict("filament_adaptive_volumetric_speed", "1");
    config.set_deserialize_strict("nozzle_temperature_initial_layer", "205");
    config.set_deserialize_strict("timelapse_type", "1");
    config.set("scan_first_layer", true);
    config.set("enable_wrapping_detection", true);
    config.set_deserialize_strict("enable_pressure_advance", "1");
    config.set_deserialize_strict("adaptive_pressure_advance", "1");
    config.set("before_layer_change_gcode", "UNSAFE_BEFORE_LAYER");
    config.set("layer_change_gcode", "UNSAFE_LAYER_CHANGE");
    config.set("time_lapse_gcode", "TIMELAPSE_TAKE_FRAME");
    config.set("change_extrusion_role_gcode", "SET_PRESSURE_ADVANCE ADVANCE=0.04");
    config.set_deserialize_strict("fan_kickstart", "0.5");
    config.set("auxiliary_fan", true);
    config.set_deserialize_strict("activate_air_filtration", "1");
    config.set_deserialize_strict("activate_chamber_temp_control", "1");

    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);

    constexpr double flow_ratio = 0.93;
    constexpr double continuous_solid_fill_ratio = 1.019;
    constexpr double target_layer_volume = 24.0 * 24.0 * 0.2;
    const double filament_area = PI * 1.75 * 1.75 * 0.25;
    bool inside = false;
    size_t sections = 0;
    double section_commanded_volume = 0.0;
    std::istringstream stream(gcode);
    std::string line;
    while (std::getline(stream, line)) {
        if (line == ";_CONTINUOUS_FERMAT_BEGIN") {
            REQUIRE_FALSE(inside);
            inside = true;
            section_commanded_volume = 0.0;
        } else if (line == ";_CONTINUOUS_FERMAT_END") {
            REQUIRE(inside);
            CHECK(section_commanded_volume / target_layer_volume ==
                  Catch::Approx(flow_ratio * continuous_solid_fill_ratio).margin(0.005));
            inside = false;
            ++sections;
        } else if (inside && line.rfind("G1", 0) == 0) {
            const auto words = parse_numeric_words(line);
            if (const auto e = words.find('E'); e != words.end())
                section_commanded_volume += e->second * filament_area;
        }
    }
    CHECK_FALSE(inside);
    CHECK(sections == 3);
    CHECK(gcode.find("\nUNSAFE_BEFORE_LAYER\n") == std::string::npos);
    CHECK(gcode.find("\nUNSAFE_LAYER_CHANGE\n") == std::string::npos);
    CHECK(gcode.find("\nTIMELAPSE_TAKE_FRAME\n") == std::string::npos);
    CHECK(gcode.find("\nSET_PRESSURE_ADVANCE ADVANCE=0.04\n") == std::string::npos);
    CHECK(gcode.find("EXCLUDE_OBJECT_START") == std::string::npos);
}

TEST_CASE("Continuous slicing export rejects a protected bead outside a nonrectangular bed", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set_key_value(
        "printable_area",
        new ConfigOptionPoints({Vec2d(0.0, 0.0), Vec2d(200.0, 0.0), Vec2d(0.0, 200.0)}));

    Print print;
    Model model;
    apply_continuous_prism(print, model, config, Vec3d(170.0, 170.0, 0.0));
    REQUIRE(print.validate().string.empty());
    REQUIRE(continuous_export_throws_slicing_error(print));
}

TEST_CASE("Continuous slicing export rejects emitted Z above printable height", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set("printable_height", 0.8);
    config.set("z_offset", 0.3);

    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    REQUIRE(continuous_export_throws_slicing_error(print));
}

TEST_CASE("Continuous slicing export respects tool 0 printable height", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set("printable_height", 10.0);
    config.set_deserialize_strict("extruder_printable_height", "0.8");
    config.set("z_offset", 0.3);

    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    REQUIRE(continuous_export_throws_slicing_error(print));
}

TEST_CASE("Continuous slicing first approach cannot Z-hop above the height limit", "[PrintGCode][continuous]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    config.set("printable_height", 0.7);
    config.set_deserialize_strict("z_hop", "5");

    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    const std::string gcode = Test::gcode(print);

    double max_motion_z = 0.0;
    std::istringstream stream(gcode);
    std::string line;
    while (std::getline(stream, line)) {
        const std::string code = line.substr(0, line.find(';'));
        if (code.rfind("G0 ", 0) != 0 && code.rfind("G1 ", 0) != 0)
            continue;
        const std::map<char, double> words = parse_numeric_words(code);
        if (const auto z = words.find('Z'); z != words.end())
            max_motion_z = std::max(max_motion_z, z->second);
    }
    CHECK(max_motion_z <= 0.7 + 1e-9);
}

TEST_CASE("Continuous slicing reports finite serialized layer material deviations", "[PrintGCode][continuous][advisory]")
{
    DynamicPrintConfig config = strict_continuous_export_config();
    Print print;
    Model model;
    apply_continuous_prism(print, model, config);
    REQUIRE(print.validate().string.empty());
    print.process();

    LayerRegion *region = print.objects().front()->layers().front()->regions().front();
    REQUIRE(region->perimeters.entities.size() == 1);
    auto *collection = dynamic_cast<ExtrusionEntityCollection *>(region->perimeters.entities.front());
    REQUIRE(collection != nullptr);
    REQUIRE(collection->entities.size() == 1);
    auto *path = dynamic_cast<ExtrusionPath *>(collection->entities.front());
    REQUIRE(path != nullptr);
    REQUIRE(path->is_continuous_fermat());
    std::fill(
        path->continuous_fermat_extrusion_multipliers.begin(),
        path->continuous_fermat_extrusion_multipliers.end(),
        0.85f);

    const std::string gcode = Test::gcode(print);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_VALIDATION_WARNING serialized_material_ratio=") != std::string::npos);
    REQUIRE(gcode.find(";_CONTINUOUS_FERMAT_END") != std::string::npos);
}

SCENARIO( "PrintGCode basic functionality", "[PrintGCode][.]") {
    GIVEN("A default configuration and a print test object") {
        WHEN("the output is executed with no support material") {
            Slic3r::Print print;
            Slic3r::Model model;
            Slic3r::Test::init_print({TestMesh::cube_20x20x20}, print, model, {
                { "layer_height",					0.2 },
                { "first_layer_height",				0.2 },
                { "first_layer_extrusion_width",	0 },
                { "gcode_comments",					true },
                { "start_gcode",					"" }
                });
            std::string gcode = Slic3r::Test::gcode(print);
            THEN("Some text output is generated.") {
                REQUIRE(gcode.size() > 0);
            }
            THEN("Exported text contains slic3r version") {
                REQUIRE(gcode.find(SLIC3R_VERSION) != std::string::npos);
            }
            //THEN("Exported text contains git commit id") {
            //    REQUIRE(gcode.find("; Git Commit") != std::string::npos);
            //    REQUIRE(gcode.find(SLIC3R_BUILD_ID) != std::string::npos);
            //}
            THEN("Exported text contains extrusion statistics.") {
                REQUIRE(gcode.find("; external perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; solid infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; top infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; support material extrusion width") == std::string::npos);
                REQUIRE(gcode.find("; first layer extrusion width") == std::string::npos);
            }
            THEN("Exported text does not contain cooling markers (they were consumed)") {
                REQUIRE(gcode.find(";_EXTRUDE_SET_SPEED") == std::string::npos);
            }

            THEN("GCode preamble is emitted.") {
                REQUIRE(gcode.find("G21 ; set units to millimeters") != std::string::npos);
            }

            THEN("Config options emitted for print config, default region config, default object config") {
                REQUIRE(gcode.find("; first_layer_temperature") != std::string::npos);
                REQUIRE(gcode.find("; layer_height") != std::string::npos);
                REQUIRE(gcode.find("; fill_density") != std::string::npos);
            }
            THEN("Infill is emitted.") {
                boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, infill_regex));
            }
            THEN("Perimeters are emitted.") {
				boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, perimeters_regex));
            }
            THEN("Skirt is emitted.") {
                boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, skirt_regex));
            }
            THEN("final Z height is 20mm") {
                double final_z = 0.0;
                GCodeReader reader;
                reader.apply_config(print.config());
                reader.parse_buffer(gcode, [&final_z] (GCodeReader& self, const GCodeReader::GCodeLine& line) {
                    final_z = std::max<double>(final_z, static_cast<double>(self.z())); // record the highest Z point we reach
                });
                REQUIRE(final_z == Catch::Approx(20.));
            }
        }
        WHEN("output is executed with complete objects and two differently-sized meshes") {
            Slic3r::Print print;
            Slic3r::Model model;
            Slic3r::Test::init_print({TestMesh::cube_20x20x20,TestMesh::cube_20x20x20}, print, model, {
                { "first_layer_extrusion_width",    0 },
                { "first_layer_height",             0.3 },
                { "layer_height",                   0.2 },
                { "support_material",               false },
                { "raft_layers",                    0 },
                { "complete_objects",               true },
                { "gcode_comments",                 true },
                { "between_objects_gcode",          "; between-object-gcode" }
                });
            std::string gcode = Slic3r::Test::gcode(print);
            THEN("Some text output is generated.") {
                REQUIRE(gcode.size() > 0);
            }
            THEN("Infill is emitted.") {
                boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, infill_regex));
            }
            THEN("Perimeters are emitted.") {
                boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, perimeters_regex));
            }
            THEN("Skirt is emitted.") {
                boost::smatch has_match;
                REQUIRE(boost::regex_search(gcode, has_match, skirt_regex));
            }
            THEN("Between-object-gcode is emitted.") {
                REQUIRE(gcode.find("; between-object-gcode") != std::string::npos);
            }
            THEN("final Z height is 20.1mm") {
                double final_z = 0.0;
                GCodeReader reader;
                reader.apply_config(print.config());
                reader.parse_buffer(gcode, [&final_z] (GCodeReader& self, const GCodeReader::GCodeLine& line) {
                    final_z = std::max(final_z, static_cast<double>(self.z())); // record the highest Z point we reach
                });
                REQUIRE(final_z == Catch::Approx(20.1));
            }
            THEN("Z height resets on object change") {
                double final_z = 0.0;
                bool reset = false;
                GCodeReader reader;
                reader.apply_config(print.config());
                reader.parse_buffer(gcode, [&final_z, &reset] (GCodeReader& self, const GCodeReader::GCodeLine& line) {
                    if (final_z > 0 && std::abs(self.z() - 0.3) < 0.01 ) { // saw higher Z before this, now it's lower
                        reset = true;
                    } else {
                        final_z = std::max(final_z, static_cast<double>(self.z())); // record the highest Z point we reach
                    }
                });
                REQUIRE(reset == true);
            }
            THEN("Shorter object is printed before taller object.") {
                double final_z = 0.0;
                bool reset = false;
                GCodeReader reader;
                reader.apply_config(print.config());
                reader.parse_buffer(gcode, [&final_z, &reset] (GCodeReader& self, const GCodeReader::GCodeLine& line) {
                    if (final_z > 0 && std::abs(self.z() - 0.3) < 0.01 ) { 
                        reset = (final_z > 20.0);
                    } else {
                        final_z = std::max(final_z, static_cast<double>(self.z())); // record the highest Z point we reach
                    }
                });
                REQUIRE(reset == true);
            }
        }
        WHEN("the output is executed with support material") {
            std::string gcode = ::Test::slice({TestMesh::cube_20x20x20}, {
                { "first_layer_extrusion_width",    0 },
                { "support_material",               true },
                { "raft_layers",                    3 },
                { "gcode_comments",                 true }
                });
            THEN("Some text output is generated.") {
                REQUIRE(gcode.size() > 0);
            }
            THEN("Exported text contains extrusion statistics.") {
                REQUIRE(gcode.find("; external perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; solid infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; top infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; support material extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; first layer extrusion width") == std::string::npos);
            }
            THEN("Raft is emitted.") {
                REQUIRE(gcode.find("; raft") != std::string::npos);
            }
        }
        WHEN("the output is executed with a separate first layer extrusion width") {
			std::string gcode = ::Test::slice({ TestMesh::cube_20x20x20 }, {
                { "first_layer_extrusion_width", "0.5" }
                });
            THEN("Some text output is generated.") {
                REQUIRE(gcode.size() > 0);
            }
            THEN("Exported text contains extrusion statistics.") {
                REQUIRE(gcode.find("; external perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; perimeters extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; solid infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; top infill extrusion width") != std::string::npos);
                REQUIRE(gcode.find("; support material extrusion width") == std::string::npos);
                REQUIRE(gcode.find("; first layer extrusion width") != std::string::npos);
            }
        }
        WHEN("Cooling is enabled and the fan is disabled.") {
			std::string gcode = ::Test::slice({ TestMesh::cube_20x20x20 }, {
				{ "cooling",                    true },
                { "disable_fan_first_layers",   5 }
                });
            THEN("GCode to disable fan is emitted."){
                REQUIRE(gcode.find("M107") != std::string::npos);
            }
        }
        WHEN("end_gcode exists with layer_num and layer_z") {
			std::string gcode = ::Test::slice({ TestMesh::cube_20x20x20 }, {
				{ "end_gcode",              "; Layer_num [layer_num]\n; Layer_z [layer_z]" },
                { "layer_height",           0.1 },
                { "first_layer_height",     0.1 }
                });
            THEN("layer_num and layer_z are processed in the end gcode") {
                REQUIRE(gcode.find("; Layer_num 199") != std::string::npos);
                REQUIRE(gcode.find("; Layer_z 20") != std::string::npos);
            }
        }
        WHEN("current_extruder exists in start_gcode") {
            {
				std::string gcode = ::Test::slice({ TestMesh::cube_20x20x20 }, {
					{ "start_gcode", "; Extruder [current_extruder]" }
                });
                THEN("current_extruder is processed in the start gcode and set for first extruder") {
                    REQUIRE(gcode.find("; Extruder 0") != std::string::npos);
                }
            }
			{
                DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
                config.set_num_extruders(4);
                config.set_deserialize_strict({
                    { "start_gcode",                    "; Extruder [current_extruder]" },
                    { "infill_extruder",                2 },
                    { "solid_infill_extruder",          2 },
                    { "perimeter_extruder",             2 },
                    { "support_material_extruder",      2 },
                    { "support_material_interface_extruder", 2 }
                });
                std::string gcode = Slic3r::Test::slice({TestMesh::cube_20x20x20}, config);
                THEN("current_extruder is processed in the start gcode and set for second extruder") {
                    REQUIRE(gcode.find("; Extruder 1") != std::string::npos);
                }
            }
        }

        WHEN("layer_num represents the layer's index from z=0") {
			std::string gcode = ::Test::slice({ TestMesh::cube_20x20x20, TestMesh::cube_20x20x20 }, {
				{ "complete_objects",               true },
                { "gcode_comments",                 true },
                { "layer_gcode",                    ";Layer:[layer_num] ([layer_z] mm)" },
                { "layer_height",                   0.1 },
                { "first_layer_height",             0.1 }
                });
			// End of the 1st object.
            std::string token = ";Layer:199 ";
			size_t pos = gcode.find(token);
			THEN("First and second object last layer is emitted") {
				// First object
				REQUIRE(pos != std::string::npos);
				pos += token.size();
				REQUIRE(pos < gcode.size());
				double z = 0;
				REQUIRE((sscanf(gcode.data() + pos, "(%lf mm)", &z) == 1));
				REQUIRE(z == Catch::Approx(20.));
				// Second object
				pos = gcode.find(";Layer:399 ", pos);
				REQUIRE(pos != std::string::npos);
				pos += token.size();
				REQUIRE(pos < gcode.size());
				REQUIRE((sscanf(gcode.data() + pos, "(%lf mm)", &z) == 1));
				REQUIRE(z == Catch::Approx(20.));
			}
        }
    }
}
