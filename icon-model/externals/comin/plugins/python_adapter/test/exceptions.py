import comin

try:
    comin.callback_get_ep_name(666)
    assert False, "No Exception was raised"
except comin.ComInError as exc:
    comin.print_info(f"Check successful: {exc}")


@comin.EP_SECONDARY_CONSTRUCTOR
def sec_ctr():
    # register_callback after primary constructor
    try:

        @comin.EP_ATM_TIMELOOP_START
        def foo():
            print("bar")

        assert False, "No Exception was raised"
    except comin.ComInError as exc:
        comin.print_info(f"Check successful: {exc}")


@comin.EP_ATM_TIMELOOP_BEFORE
def time_loop_start():
    try:
        metadata = comin.metadata(("foobar", 12))
        comin.print_info(dict(metadata))
        assert False, "No Exception was raised"
    except comin.ComInError as exc:
        comin.print_info(f"Check successful: {exc}")
