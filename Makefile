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
DOCKER ?= docker
DOCKER_QUIET ?= --quiet
ACT ?= act
ACT_QUIET ?= --quiet
ACT_EVENT ?= pull_request
ACT_WORKFLOW ?= .github/workflows/cmake.yml
ACT_PLATFORM ?= catthehacker/ubuntu:act-24.04
ACT_REPOSITORY_OWNER ?= tiagovla
ACT_REPOSITORY ?= $(ACT_REPOSITORY_OWNER)/oiseau
ACT_REV ?= $(shell git rev-parse HEAD)
ACT_LOCAL_REPOSITORY ?= $(ACT_REPOSITORY)@$(ACT_REV)=$(CURDIR)
ACT_IMAGE_CONTEXT ?= .github/docker/oiseau-baseline-builder
ACT_IMAGE ?= ghcr.io/$(ACT_REPOSITORY_OWNER)/oiseau-baseline-builder:latest
ACT_ARGS ?=

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

.PHONY: all configure build test precommit cppcheck website serve tidy act act-image clean

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

act-image:
	@$(DOCKER) build $(DOCKER_QUIET) --tag "$(ACT_IMAGE)" "$(ACT_IMAGE_CONTEXT)"

act: act-image
	@$(DOCKER) pull $(DOCKER_QUIET) "$(ACT_PLATFORM)"
	@$(ACT) $(ACT_EVENT) -W "$(ACT_WORKFLOW)" -P "ubuntu-26.04=$(ACT_PLATFORM)" --no-cache-server --pull=false $(ACT_QUIET) --var "OISEAU_IMAGE_OWNER=$(ACT_REPOSITORY_OWNER)" --var "OISEAU_REPOSITORY=$(ACT_REPOSITORY)" --local-repository "$(ACT_LOCAL_REPOSITORY)" $(ACT_ARGS)

clean:
	@echo "Removing generated build files..."
	@$(CMAKE) --build "$(BUILD_DIR)" --target clean
