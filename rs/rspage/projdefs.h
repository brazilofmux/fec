// Copyright (C) 1993-1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

typedef unsigned int       UINT32;
typedef          long long INT64;
typedef unsigned long long UINT64;

typedef unsigned long  ULONG;
typedef unsigned short USHORT;
typedef unsigned char  UBYTE;

#ifndef FALSE
    typedef int     BOOL;
#   define  FALSE   0
#   define  TRUE    (!FALSE)
#endif

#if defined(WIN32)
#define ROTR32(x,n) _lrotr(x,n)
#else
#define ONE32   0xFFFFFFFFUL
#define T32(x)  ((x) & ONE32)
#define ROTR32(v, n)   (T32((v) >> ((n)&0x1F)) | \
			T32((v) << (32 - ((n)&0x1F))))
#endif
