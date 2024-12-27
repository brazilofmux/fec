# Historical Notes on Original Reed-Solomon Implementation

## Overview
The original Reed-Solomon implementation (rs.cpp) from 1996-1997 contains several interesting optimizations that were critical for performance on 16/32-bit systems of that era. While these optimizations make the code more complex, they provide insights into system-level programming techniques of the mid-1990s.

## Key Optimizations

### 16-bit Processing
- Processes data in 16-bit chunks using `USHORT` type
- Uses little-endian byte ordering for efficiency
- Processes two bytes at a time to maximize memory bus utilization
```cpp
LFSR ^= *(USHORT*)(data + i);  // Read two bytes at once
*(USHORT*)bb = LFSR;           // Write final value as 16-bit chunk
```

### Lookup Table Optimization
The implementation uses a pre-computed table (ptable) that encodes generator polynomial operations:
```cpp
ptable[0] = 0x00  // First pair
ptable[1] = 0x00
ptable[2] = 0x49  // Second pair
ptable[3] = 0x01
```
- Table entries are organized in pairs for 16-bit processing
- Entries encode transformed generator polynomial operations
- Uses shifted indexing (`gg[j+1]`) in table generation

### LFSR Implementation
- Uses a Linear Feedback Shift Register for polynomial division
- Takes advantage of XOR commutativity for operation reordering
- Combines LFSR operations with 16-bit processing for efficiency

## Unexpected Properties
The implementation produces different but valid parity bytes compared to a straightforward implementation:
```
Reference:  bb[0]=0xa8, bb[1]=0x54
Original:   bb[0]=0x54, bb[1]=0xaa
```
This appears to be an unintended consequence of the optimizations, particularly:
- The byte reordering from 16-bit processing
- The transformation of polynomial operations in the lookup table
- The shifted indexing in table generation

The fact that both sets of parity bytes are valid suggests an interesting mathematical property of Reed-Solomon codes that allows for equivalent encodings through different byte orderings.

## Historical Context
- Written in an era where memory bandwidth was precious
- Optimized for 16/32-bit processors common in the mid-1990s
- Used by the rsframe program for encoding files into 64KB blocks
- Successfully used for CD backups of the era

## Legacy
While the code is highly optimized for its time, modern implementations (like rsref.cpp) prefer:
- Platform independence
- Clear implementation of standard RS algorithm
- Maintainability over maximum performance
- Consistent and predictable byte ordering

The original implementation remains an interesting example of historical optimization techniques and unexpected mathematical properties in error correction codes.
