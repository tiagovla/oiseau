.PHONY: all configure build test website serve tidy clean

all: build test

configure:
	@echo "Configuring project with CMake..."
	@cmake \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
		-DFETCHCONTENT_BASE_DIR=.cache/cmake \
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCMAKE_INSTALL_PREFIX=./opt \
		-DOISEAU_BUILD_BENCHMARK=ON \
		-DOISEAU_BUILD_DOCS=ON \
		-DOISEAU_BUILD_COVERAGE=ON \
		-DOISEAU_BUILD_TEST=ON \
		-DOISEAU_BUILD_SHARED=ON \
		-S . -B build -G Ninja

build:
	@if [ ! -f build/CMakeCache.txt ]; then $(MAKE) configure; fi
	@echo "Building project..."
	@cmake --build build --parallel

test: build
	@echo "Running tests..."
	@ctest --test-dir=build --parallel --progress --stop-on-failure --output-on-failure

website:
	@echo "Generating website..."
	cmake --build build --target website

serve:
	sleep 1 && xdg-open http://localhost:8080 &
	python -m http.server -d build/web 8080

tidy:
	find src -name '*.cpp' | xargs -r -n 1 -P $(shell nproc) clang-tidy -p build --header-filter='^src/oiseau' 2> /dev/null

clean:
	@echo "Removing generated build and install files..."
	@cmake -E rm -rf build opt
