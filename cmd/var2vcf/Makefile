.PHONY: all build test clean install help

# Binary name
BINARY_NAME=var2vcf
VERSION=1.0.0

# Build flags
LDFLAGS=-ldflags "-s -w -X main.Version=$(VERSION)"

## all: Build the binary
all: build

## build: Build the binary
build:
	@echo "Building $(BINARY_NAME)..."
	@go build $(LDFLAGS) -o $(BINARY_NAME) .
	@echo "✓ Build complete: ./$(BINARY_NAME)"

## test: Run tests
test:
	@echo "Running tests..."
	@go test -v -cover ./...

## bench: Run benchmarks
bench:
	@echo "Running benchmarks..."
	@go test -bench=. -benchmem ./...

## clean: Remove build artifacts
clean:
	@echo "Cleaning..."
	@rm -f $(BINARY_NAME)
	@go clean
	@echo "✓ Clean complete"

## install: Install the binary to $GOPATH/bin
install: build
	@echo "Installing to $$GOPATH/bin..."
	@go install $(LDFLAGS)
	@echo "✓ Installed $(BINARY_NAME)"

## validate: Validate against Perl version (requires test data)
validate: build
	@echo "Running validation against Perl version..."
	@if [ ! -f test_input.tsv ]; then \
		echo "Error: test_input.tsv not found"; \
		exit 1; \
	fi
	@echo "Testing Go version..."
	@cat test_input.tsv | ./$(BINARY_NAME) -N sample -f 0.01 > go_output.vcf
	@if command -v perl &> /dev/null && [ -f var2vcf_valid.pl ]; then \
		echo "Testing Perl version..."; \
		cat test_input.tsv | perl var2vcf_valid.pl -N sample -f 0.01 > perl_output.vcf; \
		echo "Comparing outputs (excluding headers)..."; \
		grep -v "^##source" go_output.vcf > go_variants.txt; \
		grep -v "^##source" perl_output.vcf > perl_variants.txt; \
		diff -u perl_variants.txt go_variants.txt && echo "✓ Outputs match!" || echo "✗ Outputs differ"; \
		rm -f perl_variants.txt go_variants.txt; \
	else \
		echo "Perl or var2vcf_valid.pl not found, skipping comparison"; \
	fi
	@rm -f go_output.vcf perl_output.vcf

## format: Format Go code
format:
	@echo "Formatting code..."
	@go fmt ./...
	@echo "✓ Format complete"

## lint: Run linters (requires golangci-lint)
lint:
	@echo "Running linters..."
	@if command -v golangci-lint &> /dev/null; then \
		golangci-lint run; \
	else \
		echo "golangci-lint not found, skipping"; \
	fi

## help: Show this help message
help:
	@echo "VarDict var2vcf Go Implementation"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk '/^##/ {sub(/^## /, "", $$0); printf "  %-20s %s\n", $$1, substr($$0, index($$0, ":")+2)}' $(MAKEFILE_LIST)

.DEFAULT_GOAL := help
