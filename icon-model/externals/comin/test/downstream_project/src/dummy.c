#include <stdio.h>

#include <comin.h>

static void dummy_plugin_secondary_constructor() {
  printf("DUMMY: secondary_constructor");
}

void comin_main() {
  printf("DUMMY: comin_main");
  comin_callback_register(EP_SECONDARY_CONSTRUCTOR,
                          dummy_plugin_secondary_constructor);
}
