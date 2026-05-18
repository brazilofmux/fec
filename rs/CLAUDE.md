# Reed-Solomon Library (RS) Guidelines

## Implementations
- **rskarn.cpp**: Slightly modified version of Phil Karn's 1996 implementation
- **rs.cpp**: Original adaption from 25+ years ago with extreme optimizations
- **rsref.cpp**: LLM-created reference implementation with improved readability over Karn's code
- **librs/**: Factory-based implementation offering both readability and optimal performance

All implementations produce identical results and are verified against each other using the compare_encoders tool.

## Benchmarks
Per-platform `rstest -b` results (scalar / NEON / AVX2 / AVX-512 µs/decode, plus Auto-phase breakdown) are recorded in `BENCHMARKS.md`. Append a new dated section there when you re-run on a new box rather than editing comments in headers.

## Build Commands
- Build library: `cd librs && make`
- Build tests: `cd rstest && make`
- Run benchmarks: `cd rstest && ./rstest -b`
- Run comparison test: `cd rstest && ./rstest -c 6` (tests tt=6)
- Run encoding test: `cd rstest && ./rstest -e 10 -n 1000` (1000 codewords with tt=10)
- Run decoding test: `cd rstest && ./rstest -d 10 -n 1000` (decode 1000 codewords with tt=10)
- Compare encoders: `cd rstest && ./compare_encoders`
- Build/run frame code: `cd rsframe && make && ./rsframe`

## Code Style Guidelines
- **Naming**: Classes use PascalCase, variables use snake_case, constants use UPPER_CASE
- **Formatting**: 4-space indentation, braces on same line for functions
- **Headers**: Use include guards with uppercase project prefix
- **Error Handling**: Return integer error codes, check return values
- **Types**: Use custom type aliases for domain-specific types (e.g., `typedef unsigned char GF`)
- **Organization**: Base classes with specialized implementations
- **Comments**: Document mathematical algorithms and time complexity
- **Testing**: All encoders/decoders must produce identical results for identical inputs

Always maintain backward compatibility with existing interfaces.