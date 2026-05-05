import comin


@comin.EP_SECONDARY_CONSTRUCTOR
def sec_ctr():
    raise RuntimeError("FOOO")
