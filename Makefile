.DEFAULT_GOAL := all
SHELL := /bin/sh

# Build configuration
BUILD_DIR ?= build
BUILD_TYPE ?= Debug
CMAKE_GENERATOR ?= Ninja
JOBS ?= $(shell nproc)

# Tool commands
CMAKE ?= cmake
CTEST ?= ctest
PRE_COMMIT ?= pre-commit
CPPCHECK ?= cppcheck
CLANG_TIDY ?= clang-tidy
PYTHON ?= python3
OPEN_BROWSER ?= xdg-open

# Website serving
HTTP_PORT ?= 8080
BROWSER_OPEN_DELAY ?= 1

# Project options
OISEAU_BUILD_DEMO ?= ON
OISEAU_BUILD_TEST ?= ON
OISEAU_BUILD_DOCS ?= ON
OISEAU_BUILD_BENCHMARK ?= ON
OISEAU_BUILD_COVERAGE ?= OFF
OISEAU_COVERAGE_AUTO_INSTALL_GCOVR ?= ON
OISEAU_BUILD_SHARED ?= ON
OISEAU_PRECOMPILE_HEADERS ?= ON

# Tool options
CTEST_FLAGS ?= --parallel --progress --stop-on-failure --output-on-failure
CTEST_COLOR ?= 1
CPPCHECK_FLAGS ?= \
	--enable=all \
	--std=c++20 \
	--suppress=checkersReport \
	--suppress=missingInclude \
	--suppress=unmatchedSuppression \
	--suppress=missingIncludeSystem \
	--suppress=unusedFunction \
	--suppress=returnByReference \
	--error-exitcode=2 \
	--check-level=exhaustive
CLANG_TIDY_FLAGS ?= --header-filter='^src/oiseau'

CMAKE_ARGS ?= \
	-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
	-DOISEAU_BUILD_DEMO=$(OISEAU_BUILD_DEMO) \
	-DOISEAU_BUILD_TEST=$(OISEAU_BUILD_TEST) \
	-DOISEAU_BUILD_DOCS=$(OISEAU_BUILD_DOCS) \
	-DOISEAU_BUILD_BENCHMARK=$(OISEAU_BUILD_BENCHMARK) \
	-DOISEAU_BUILD_COVERAGE=$(OISEAU_BUILD_COVERAGE) \
	-DOISEAU_COVERAGE_AUTO_INSTALL_GCOVR=$(OISEAU_COVERAGE_AUTO_INSTALL_GCOVR) \
	-DOISEAU_BUILD_SHARED=$(OISEAU_BUILD_SHARED) \
	-DOISEAU_PRECOMPILE_HEADERS=$(OISEAU_PRECOMPILE_HEADERS)
CMAKE_EXTRA_ARGS ?=

.PHONY: all configure build test precommit cppcheck website serve tidy clean

all: build test precommit cppcheck tidy

configure:
	@echo "Configuring project with CMake..."
	@$(CMAKE) -S . -B "$(BUILD_DIR)" -G "$(CMAKE_GENERATOR)" $(CMAKE_ARGS) $(CMAKE_EXTRA_ARGS)

build: configure
	@echo "Building project..."
	@$(CMAKE) --build "$(BUILD_DIR)" --parallel "$(JOBS)"

test: build
	@echo "Running tests..."
	@NO_COLOR= CLICOLOR_FORCE=$(CTEST_COLOR) $(CTEST) --test-dir "$(BUILD_DIR)" $(CTEST_FLAGS)

precommit:
	@echo "Running pre-commit hooks..."
	@$(PRE_COMMIT) run --all-files

cppcheck:
	@echo "Running cppcheck..."
	@$(CPPCHECK) $(CPPCHECK_FLAGS) src -I src

website: build
	@echo "Generating website..."
	@$(CMAKE) --build "$(BUILD_DIR)" --target website --parallel "$(JOBS)"

serve: website
	@sleep "$(BROWSER_OPEN_DELAY)" && $(OPEN_BROWSER) "http://localhost:$(HTTP_PORT)" &
	@$(PYTHON) -m http.server --directory "$(BUILD_DIR)/web" "$(HTTP_PORT)"

tidy: configure
	@find src -name '*.cpp' -print0 | xargs -0 -r -n 1 -P "$(JOBS)" $(CLANG_TIDY) -p "$(BUILD_DIR)" $(CLANG_TIDY_FLAGS) 2>/dev/null

clean:
	@echo "Removing generated build files..."
	@$(CMAKE) --build "$(BUILD_DIR)" --target clean
