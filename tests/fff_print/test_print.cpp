#include <catch2/catch_all.hpp>

#include "libslic3r/libslic3r.h"
#include "libslic3r/Print.hpp"
#include "libslic3r/Layer.hpp"

#include "test_data.hpp"

#include <functional>
#include <limits>

using namespace Slic3r;
using namespace Slic3r::Test;

namespace {

DynamicPrintConfig strict_continuous_validation_config()
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"spiral_mode", true},
        {"spiral_hybrid_non_crossing", true},
        {"gcode_flavor", "marlin2"},
        {"print_sequence", "by layer"},
        {"use_relative_e_distances", true},
        {"enable_power_loss_recovery", "disable"},
        {"exclude_object", false},
        {"nozzle_temperature_initial_layer", "200"},
        {"nozzle_temperature", "200"},
        {"filament_flow_ratio", "1"},
        {"filament_max_volumetric_speed", "12"},
        {"filament_adaptive_volumetric_speed", "0"},
        {"max_volumetric_extrusion_rate_slope", 0.0},
        {"enable_pressure_advance", "0"},
        {"adaptive_pressure_advance", "0"},
        {"emit_machine_limits_to_gcode", false},
        {"gcode_add_line_number", false},
        {"file_start_gcode", ""},
        {"machine_start_gcode", ""},
        {"filament_start_gcode", ""},
        {"machine_end_gcode", ""},
        {"filament_end_gcode", ""},
        {"enable_support", false},
        {"enforce_support_layers", 0},
        {"raft_layers", 0},
        {"skirt_loops", 0},
        {"skirt_height", 0},
        {"draft_shield", "disabled"},
        {"brim_type", "no_brim"},
        {"brim_width", 0.0},
        {"enable_prime_tower", false},
        {"timelapse_type", "0"},
        {"enable_wrapping_detection", false},
    });
    return config;
}

void apply_cubes_without_arranging(
    Print &print,
    Model &model,
    const DynamicPrintConfig &config,
    const size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        ModelObject *object = model.add_object();
        object->name = "continuous-test-cube.stl";
        object->add_volume(Test::mesh(TestMesh::cube_20x20x20));
        object->add_instance();
        object->ensure_on_bed();
        print.auto_assign_extruders(object);
    }
    print.apply(model, config);
}

StringObjectException validate_cube(const DynamicPrintConfig &config)
{
    Print print;
    Model model;
    apply_cubes_without_arranging(print, model, config, 1);
    return print.validate();
}

} // namespace

SCENARIO("Print: strict continuous slicing validation is fail-closed", "[Print][validate][continuous]")
{
    GIVEN("one object, one instance, one material, and no auxiliary extrusion")
    {
        DynamicPrintConfig config = strict_continuous_validation_config();

        THEN("the strict configuration is accepted")
        {
            CHECK(validate_cube(config).string.empty());
        }

        struct InvalidFeature {
            const char *name;
            std::function<void(DynamicPrintConfig &)> enable;
        };
        const std::vector<InvalidFeature> invalid_features{
            {"support material", [](DynamicPrintConfig &cfg) { cfg.set("enable_support", true); }},
            {"a raft", [](DynamicPrintConfig &cfg) { cfg.set("raft_layers", 1); }},
            {"a skirt", [](DynamicPrintConfig &cfg) {
                 cfg.set("skirt_loops", 1);
                 cfg.set("skirt_height", 1);
             }},
            {"a draft shield", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("draft_shield", "enabled"); }},
            {"a brim", [](DynamicPrintConfig &cfg) {
                 cfg.set_deserialize_strict("brim_type", "outer_only");
                 cfg.set("brim_width", 2.0);
             }},
            {"smooth timelapse", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("timelapse_type", "1"); }},
            {"first-layer scanning", [](DynamicPrintConfig &cfg) { cfg.set("scan_first_layer", true); }},
            {"clumping detection", [](DynamicPrintConfig &cfg) { cfg.set("enable_wrapping_detection", true); }},
            {"ironing", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("ironing_type", "top"); }},
            {"before-layer executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("before_layer_change_gcode", "G1 X10"); }},
            {"layer-change executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("layer_change_gcode", "G92 E0"); }},
            {"timelapse executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("time_lapse_gcode", "G1 X10"); }},
            {"role-change executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("change_extrusion_role_gcode", "G1 E1"); }},
            {"bare-CR hook bypass", [](DynamicPrintConfig &cfg) { cfg.set("layer_change_gcode", "; comment\rG1 X10"); }},
            {"non-Marlin G-code semantics", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("gcode_flavor", "klipper"); }},
            {"Marlin 1 G-code semantics", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("gcode_flavor", "marlin"); }},
            {"the continuous flag without spiral mode", [](DynamicPrintConfig &cfg) { cfg.set("spiral_mode", false); }},
            {"absolute extrusion distances", [](DynamicPrintConfig &cfg) { cfg.set("use_relative_e_distances", false); }},
            {"object cancellation", [](DynamicPrintConfig &cfg) { cfg.set("exclude_object", true); }},
            {"power-loss recovery", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("enable_power_loss_recovery", "enable"); }},
            {"pressure equalization", [](DynamicPrintConfig &cfg) { cfg.set("max_volumetric_extrusion_rate_slope", 1.0); }},
            {"fixed pressure advance", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("enable_pressure_advance", "1"); }},
            {"adaptive pressure advance", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("adaptive_pressure_advance", "1"); }},
            {"firmware machine-limit overrides", [](DynamicPrintConfig &cfg) { cfg.set("emit_machine_limits_to_gcode", true); }},
            {"post-export line numbering", [](DynamicPrintConfig &cfg) { cfg.set("gcode_add_line_number", true); }},
            {"file-start executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("file_start_gcode", "G28"); }},
            {"machine-start executable G-code", [](DynamicPrintConfig &cfg) { cfg.set("machine_start_gcode", "G1 X10"); }},
            {"filament-start executable G-code", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("filament_start_gcode", "G1 E5"); }},
            {"a nonzero physical tool mapping", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("physical_extruder_map", "1"); }},
            {"a nozzle temperature transition", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("nozzle_temperature_initial_layer", "210"); }},
            {"a non-unity filament flow ratio", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("filament_flow_ratio", "1.1"); }},
            {"no volumetric speed limit", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("filament_max_volumetric_speed", "0"); }},
            {"an adaptive volumetric limit", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("filament_adaptive_volumetric_speed", "1"); }},
            {"a non-positive printable height", [](DynamicPrintConfig &cfg) { cfg.set("printable_height", 0.0); }},
            {"a non-finite Z offset", [](DynamicPrintConfig &cfg) { cfg.set("z_offset", std::numeric_limits<double>::quiet_NaN()); }},
            {"a bed exclusion polygon", [](DynamicPrintConfig &cfg) {
                 cfg.set_key_value(
                     "bed_exclude_area",
                     new ConfigOptionPoints({Vec2d(5.0, 5.0), Vec2d(15.0, 5.0), Vec2d(15.0, 15.0), Vec2d(5.0, 15.0)}));
             }},
            {"a tool-specific printable polygon", [](DynamicPrintConfig &cfg) {
                 cfg.set_key_value(
                     "extruder_printable_area",
                     new ConfigOptionPointsGroups({{Vec2d(0.0, 0.0), Vec2d(100.0, 0.0), Vec2d(100.0, 100.0), Vec2d(0.0, 100.0)}}));
             }},
            {"sequential printing", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("print_sequence", "by object"); }},
            {"machine end motion", [](DynamicPrintConfig &cfg) { cfg.set("machine_end_gcode", "G28 X0"); }},
            {"fan kickstart motion splitting", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("fan_kickstart", "0.5"); }},
            {"an auxiliary fan", [](DynamicPrintConfig &cfg) { cfg.set("auxiliary_fan", true); }},
            {"air filtration", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("activate_air_filtration", "1"); }},
            {"chamber control", [](DynamicPrintConfig &cfg) { cfg.set_deserialize_strict("activate_chamber_temp_control", "1"); }},
        };

        for (const InvalidFeature &feature : invalid_features) {
            DYNAMIC_SECTION("rejects " << feature.name)
            {
                feature.enable(config);
                const StringObjectException error = validate_cube(config);
                CHECK_FALSE(error.string.empty());
                CHECK(error.opt_key == "spiral_hybrid_non_crossing");
            }
        }
    }
}

SCENARIO("Print: strict continuous slicing rejects an otherwise unused second nozzle", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    Print print;
    Model model;
    apply_cubes_without_arranging(print, model, config, 1);
    const_cast<PrintConfig &>(print.config()).nozzle_diameter.values.push_back(0.4);

    const StringObjectException error = print.validate();
    CHECK_FALSE(error.string.empty());
    CHECK(error.opt_key == "spiral_hybrid_non_crossing");
}

SCENARIO("Print: strict continuous slicing rejects Bambu-specific output", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    Print print;
    Model model;
    print.is_BBL_printer() = true;
    apply_cubes_without_arranging(print, model, config, 1);

    const StringObjectException error = print.validate();
    CHECK_FALSE(error.string.empty());
    CHECK(error.opt_key == "spiral_hybrid_non_crossing");
}

SCENARIO("Print: strict continuous slicing requires a deterministic finite plate origin", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    Print print;
    Model model;
    CHECK(print.get_plate_origin().isZero());
    print.set_plate_origin(Vec3d(0.0, std::numeric_limits<double>::quiet_NaN(), 0.0));
    apply_cubes_without_arranging(print, model, config, 1);

    const StringObjectException error = print.validate();
    CHECK_FALSE(error.string.empty());
    CHECK(error.opt_key == "spiral_hybrid_non_crossing");
}

SCENARIO("Print: strict continuous slicing requires one object and one instance", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    Print print;
    Model model;

    GIVEN("two print objects")
    {
        apply_cubes_without_arranging(print, model, config, 2);
        const StringObjectException error = print.validate();
        CHECK_FALSE(error.string.empty());
        CHECK(error.opt_key == "spiral_hybrid_non_crossing");
    }

    GIVEN("two instances of one print object")
    {
        apply_cubes_without_arranging(print, model, config, 1);
        model.objects.front()->add_instance();
        print.apply(model, config);
        const StringObjectException error = print.validate();
        CHECK_FALSE(error.string.empty());
        CHECK(error.opt_key == "spiral_hybrid_non_crossing");
    }
}

SCENARIO("Print: strict continuous slicing permits comment-only layer hooks", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    config.set("before_layer_change_gcode", "; before {layer_num}");
    config.set("layer_change_gcode", "  ; layer {layer_num}");
    config.set("time_lapse_gcode", "; camera intentionally disabled");
    config.set("change_extrusion_role_gcode", "; role {extrusion_role}");
    CHECK(validate_cube(config).string.empty());
}

SCENARIO("Print: classic spiral mode keeps its existing validation", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    config.set("spiral_hybrid_non_crossing", false);
    config.set("use_relative_e_distances", false);
    config.set("skirt_loops", 1);
    config.set("skirt_height", 1);
    config.set_deserialize_strict("brim_type", "outer_only");
    config.set("brim_width", 2.0);

    CHECK(validate_cube(config).string.empty());
}

SCENARIO("Print: toggling continuous slicing invalidates the object slice", "[Print][validate][continuous]")
{
    DynamicPrintConfig config = strict_continuous_validation_config();
    config.set("spiral_hybrid_non_crossing", false);

    Print print;
    Model model;
    apply_cubes_without_arranging(print, model, config, 1);
    print.process();
    REQUIRE(print.objects().front()->is_step_done(posSlice));

    config.set("spiral_hybrid_non_crossing", true);
    print.apply(model, config);
    CHECK_FALSE(print.objects().front()->is_step_done(posSlice));
}

SCENARIO("PrintObject: Perimeter generation", "[PrintObject][.]") {
    GIVEN("20mm cube and default config") {
        WHEN("make_perimeters() is called")  {
            Slic3r::Print print;
            Slic3r::Test::init_and_process_print({TestMesh::cube_20x20x20}, print, { { "fill_density", 0 } });
			const PrintObject &object = *print.objects().front();
			THEN("67 layers exist in the model") {
                REQUIRE(object.layers().size() == 66);
            }
            THEN("Every layer in region 0 has 1 island of perimeters") {
                for (const Layer *layer : object.layers())
                    REQUIRE(layer->regions().front()->perimeters.entities.size() == 1);
            }
            THEN("Every layer in region 0 has 3 paths in its perimeters list.") {
                for (const Layer *layer : object.layers())
                    REQUIRE(layer->regions().front()->perimeters.items_count() == 3);
            }
        }
    }
}

SCENARIO("Print: Skirt generation", "[Print][.]") {
    GIVEN("20mm cube and default config") {
        WHEN("Skirts is set to 2 loops")  {
            Slic3r::Print print;
            Slic3r::Test::init_and_process_print({TestMesh::cube_20x20x20}, print, {
            	{ "skirt_height", 	1 },
        		{ "skirt_distance", 1 },
        		{ "skirts", 		2 }
            });
            THEN("Skirt Extrusion collection has 2 loops in it") {
                REQUIRE(print.skirt().items_count() == 2);
                REQUIRE(print.skirt().flatten().entities.size() == 2);
            }
        }
    }
}

SCENARIO("Print: Changing number of solid surfaces does not cause all surfaces to become internal.", "[Print][.]") {
    GIVEN("sliced 20mm cube and config with top_solid_surfaces = 2 and bottom_solid_surfaces = 1") {
        Slic3r::DynamicPrintConfig config = Slic3r::DynamicPrintConfig::full_print_config();
		config.set_deserialize_strict({
			{ "top_solid_layers",		2 },
			{ "bottom_solid_layers",	1 },
			{ "layer_height",			0.25 }, // get a known number of layers
			{ "first_layer_height",		0.25 }
			});
        Slic3r::Print print;
        Slic3r::Model model;
        Slic3r::Test::init_print({TestMesh::cube_20x20x20}, print, model, config);
        // Precondition: Ensure that the model has 2 solid top layers (39, 38)
        // and one solid bottom layer (0).
		auto test_is_solid_infill = [&print](size_t obj_id, size_t layer_id) {
		    const Layer &layer = *(print.objects().at(obj_id)->get_layer((int)layer_id));
		    // iterate over all of the regions in the layer
		    for (const LayerRegion *region : layer.regions()) {
		        // for each region, iterate over the fill surfaces
		        for (const Surface &surface : region->fill_surfaces.surfaces)
		            CHECK(surface.is_solid());
		    }
		};
        print.process();
        test_is_solid_infill(0,  0); // should be solid
        test_is_solid_infill(0, 79); // should be solid
        test_is_solid_infill(0, 78); // should be solid
        WHEN("Model is re-sliced with top_solid_layers == 3") {
			config.set("top_solid_layers", 3);
			print.apply(model, config);
            print.process();
            THEN("Print object does not have 0 solid bottom layers.") {
                test_is_solid_infill(0, 0);
            }
            AND_THEN("Print object has 3 top solid layers") {
                test_is_solid_infill(0, 79);
                test_is_solid_infill(0, 78);
                test_is_solid_infill(0, 77);
            }
        }
    }
}

SCENARIO("Print: Brim generation", "[Print][.]") {
    GIVEN("20mm cube and default config, 1mm first layer width") {
        WHEN("Brim is set to 3mm")  {
	        Slic3r::Print print;
	        Slic3r::Test::init_and_process_print({TestMesh::cube_20x20x20}, print, {
	        	{ "first_layer_extrusion_width", 	1 },
	        	{ "brim_width", 					3 }
	        });
            THEN("Brim Extrusion collection has 3 loops in it") {
                size_t total_items = 0;
                for (const auto& pair : print.get_brimMap()) {
                    total_items += pair.second.items_count();
                }
                REQUIRE(total_items == 3);
            }
        }
        WHEN("Brim is set to 6mm")  {
	        Slic3r::Print print;
	        Slic3r::Test::init_and_process_print({TestMesh::cube_20x20x20}, print, {
	        	{ "first_layer_extrusion_width", 	1 },
	        	{ "brim_width", 					6 }
	        });
            THEN("Brim Extrusion collection has 6 loops in it") {
                size_t total_items = 0;
                for (const auto& pair : print.get_brimMap()) {
                    total_items += pair.second.items_count();
                }
                REQUIRE(total_items == 6);
            }
        }
        WHEN("Brim is set to 6mm, extrusion width 0.5mm")  {
	        Slic3r::Print print;
	        Slic3r::Test::init_and_process_print({TestMesh::cube_20x20x20}, print, {
	        	{ "first_layer_extrusion_width", 	1 },
	        	{ "brim_width", 					6 },
	        	{ "first_layer_extrusion_width", 	0.5 }
	        });
			print.process();
            THEN("Brim Extrusion collection has 12 loops in it") {
                size_t total_items = 0;
                for (const auto& pair : print.get_brimMap()) {
                    total_items += pair.second.items_count();
                }
                REQUIRE(total_items == 14);
            }
        }
    }
}
