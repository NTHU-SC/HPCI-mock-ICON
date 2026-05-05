# ICON RAGNAROK

C++ implementation of the AES Physics based on [Kokkos](https://github.com/kokkos/kokkos) Core programming model.

## Usage

The project is part of ICON and built by the ICON build system if the option `--enable-ragnarok` is turned on; it also handles the proper compilation flags needed to configure Kokkos backend (from `externals/kokkos` submodule).

## Available compile options (default options are marked in __bold__)
* _Unit-test_ - compile unit-tests
  * BUILD_TESTING=__ON__/OFF
* Standalone - compile standalone components
  * RGK_ENABLE_STANDALONE=__ON__/OFF

### Build as part of ICON

ICON configuration script will also configure ragnarok & kokkos.

### Build standalone

To build ragnarok standalone, one can use the cache files provided for levante under `ragnarok/cmake/caches/`, e.g.:

  ```
  cmake -B build -S ragnarok -C ragnarok/cmake/caches/levante_cpu.cmake`
  ```

## Testing

The unit tests support via gtest. To run the unit-tests, e.g.:

  ```
  cd build && ctest
  ```
