BUILD_DIR := build

.PHONY: all configure build test precommit cppcheck website serve tidy clean

all: build test precommit cppcheck

configure:
	@echo "Configuring project with CMake..."
	@cmake -S . -B $(BUILD_DIR) -G Ninja \
		-DCMAKE_BUILD_TYPE=Debug \
		-DOISEAU_BUILD_BENCHMARK=ON \
		-DOISEAU_BUILD_DOCS=ON \
		-DOISEAU_BUILD_COVERAGE=ON \
		-DOISEAU_BUILD_TEST=ON \
		-DOISEAU_BUILD_SHARED=ON

build:
	@echo "Building project..."
	@cmake --build $(BUILD_DIR) --parallel

test: build
	@echo "Running tests..."
	@ctest --test-dir=$(BUILD_DIR) --parallel --progress --stop-on-failure --output-on-failure

precommit:
	@echo "Running pre-commit hooks..."
	@pre-commit run --all-files

cppcheck:
	@echo "Running cppcheck..."
	@cppcheck \
		--enable=all \
		--std=c++20 \
		--suppress=missingInclude \
		--suppress=unmatchedSuppression \
		--suppress=missingIncludeSystem \
		--suppress=unusedFunction \
		--suppress=returnByReference \
		--error-exitcode=2 \
		--check-level=exhaustive \
		src \
		-I src

website:
	@echo "Generating website..."
	cmake --build $(BUILD_DIR) --target website

serve:
	sleep 1 && xdg-open http://localhost:8080 &
	python -m http.server -d $(BUILD_DIR)/web 8080

tidy:
	find src -name '*.cpp' | xargs -r -n 1 -P $(shell nproc) clang-tidy -p $(BUILD_DIR) --header-filter='^src/oiseau' 2> /dev/null

clean:
	@echo "Removing generated build and install files..."
	@cmake --build $(BUILD_DIR) --target clean
