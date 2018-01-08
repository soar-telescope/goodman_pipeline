from __future__ import absolute_import


# import all classes in core.py
from ..core import (CriticalError,
                    NightDataContainer,
                    NoMatchFound,
                    NotEnoughLinesDetected,
                    NoTargetException,
                    SpectroscopicMode)

# import of functions in core.py
from ..core import (add_wcs_keys,
                    call_cosmic_rejection,
                    classify_spectroscopic_data,
                    convert_time,
                    create_background_image,
                    dcr_cosmicray_rejection,
                    extract,
                    extract_fractional_pixel,
                    extract_optimal,
                    extract_simple,
                    fix_duplicated_keywords,
                    fractional_sum,
                    get_background_value,
                    get_best_flat,
                    get_central_wavelength,
                    get_extraction_zone,
                    get_slit_trim_section,
                    get_twilight_time,
                    identify_targets,
                    image_overscan,
                    image_trim,
                    lacosmic_cosmicray_rejection,
                    normalize_master_flat,
                    print_default_args,
                    print_progress,
                    print_spacers,
                    ra_dec_to_deg,
                    read_fits,
                    remove_background_by_median,
                    remove_conflictive_keywords,
                    search_comp_group,
                    spectroscopic_extraction,
                    trace,
                    trace_targets,
                    write_fits)


def test_critical_error():
    pass


def test_night_data_container():
    pass


def test_no_match_found():
    pass


def test_not_enough_lines_detected():
    pass


def test_no_target_exception():
    pass


def test_spectroscopic_mode():
    pass


def test_add_wcs_keys():
    pass


def test_call_cosmic_rejection():
    pass


def test_classify_spectroscopic_data():
    pass


def test_convert_time():
    pass


def test_create_background_image():
    pass


def test_dcr_cosmicray_rejection():
    pass


def test_extract():
    pass


def test_extract_fractional_pixel():
    pass


def test_extract_optimal():
    pass


def test_extract_simple():
    pass


def test_fix_duplicated_keywords():
    pass


def test_fractional_sum():
    pass


def test_get_background_value():
    pass


def test_get_best_flat():
    pass


def test_get_central_wavelength():
    pass


def test_get_extraction_zone():
    pass


def test_get_slit_trim_section():
    pass


def test_get_twilight_time():
    pass


def test_identify_targets():
    pass


def test_image_overscan():
    pass


def test_image_trim():
    pass


def test_lacosmic_cosmicray_rejection():
    pass


def test_normalize_master_flat():
    pass


def test_print_default_args():
    pass


def test_print_progress():
    pass


def test_print_spacers():
    pass


def test_ra_dec_to_deg():
    pass


def test_read_fits():
    pass


def test_remove_background_by_median():
    pass


def test_remove_conflictive_keywords():
    pass


def test_search_comp_group():
    pass


def test_spectroscopic_extraction():
    pass


def test_trace():
    pass


def test_trace_targets():
    pass


def test_write_fits():
    pass
