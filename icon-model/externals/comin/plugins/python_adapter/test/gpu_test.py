import comin
import numpy as np

glob = comin.descrdata_get_global()
if glob.has_device:
    comin.print_info(f"{glob.device_name=}")
    comin.print_info(f"{glob.device_vendor=}")
    comin.print_info(f"{glob.device_driver=}")

if glob.has_device and "NVIDIA" in glob.device_vendor.upper():
    try:
        comin.print_info("Using cupy!")
        import cupy as xp

        DEVICE_SYNC_FLAG = comin.COMIN_FLAG_DEVICE
    except ImportError as e:
        comin.print_info("Cannot import cupy, falling back to numpy")
        comin.print_info(e)
        import sys

        comin.print_info(sys.path)
        import numpy as xp

        DEVICE_SYNC_FLAG = 0
else:
    comin.print_info("No NVIDIA device found falling back to numpy")
    import numpy as xp

    DEVICE_SYNC_FLAG = 0

domain = comin.descrdata_get_domain(1)

comin.var_request_add(("device_to_host", 1), True)
comin.var_request_add(("host_to_device", 1), True)


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def sec_ctor():
    global ta, device_to_host1, device_to_host2, host_to_device1, host_to_device2
    ta = comin.var_get(
        [comin.EP_ATM_WRITE_OUTPUT_BEFORE], ("temp", 1), comin.COMIN_FLAG_READ
    )
    device_to_host1 = comin.var_get(
        [comin.EP_ATM_NUDGING_BEFORE],
        ("device_to_host", 1),
        comin.COMIN_FLAG_WRITE | DEVICE_SYNC_FLAG,
    )
    host_to_device1 = comin.var_get(
        [comin.EP_ATM_NUDGING_BEFORE],
        ("host_to_device", 1),
        comin.COMIN_FLAG_WRITE,
    )
    device_to_host2 = comin.var_get(
        [comin.EP_ATM_NUDGING_AFTER], ("device_to_host", 1), comin.COMIN_FLAG_READ
    )
    host_to_device2 = comin.var_get(
        [comin.EP_ATM_NUDGING_AFTER],
        ("host_to_device", 1),
        comin.COMIN_FLAG_READ | DEVICE_SYNC_FLAG,
    )


@comin.register_callback(comin.EP_ATM_WRITE_OUTPUT_BEFORE)
def foo():
    comin.print_info(f"{ta.__cuda_array_interface__=}")
    ta_arr = np.asarray(ta)
    if hasattr(ta_arr, "__cuda_array_interface__"):
        comin.print_info(f"{ta_arr.__cuda_array_interface__=}")
    comin.print_info(f"{type(ta_arr)=}")
    if hasattr(ta_arr, "device"):
        comin.print_info(f"{ta_arr.device=}")
    comin.print_info(f"{ta_arr.base}")
    comin.print_info("Computing mean surface temperture (on this process)")
    tas = ta_arr[:, -1, :, 0, 0]
    comin.print_info(f"{tas.mean()=}")


@comin.register_callback(comin.EP_ATM_NUDGING_BEFORE)
def set_to_42():
    device_to_host_xp = xp.asarray(device_to_host1)
    device_to_host_xp[:] = 42.0
    host_to_device_np = np.asarray(host_to_device1)
    host_to_device_np[:] = 43.0


@comin.register_callback(comin.EP_ATM_NUDGING_AFTER)
def print_element():
    device_to_host_np = np.asarray(device_to_host2)
    assert np.allclose(device_to_host_np, 42.0)
    comin.print_info("check successful for device_to_host")
    host_to_device_xp = xp.asarray(host_to_device2)
    assert xp.allclose(host_to_device_xp, 43.0)
    comin.print_info("check successful for host_to_device")
