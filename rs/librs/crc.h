// Copyright (C) 1993-1996 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

#include <stdint.h>

#define CRC32_INITIAL        0xFFFFFFFFU
#define CRC32_VALID_RESULT   0x2144DF1CU

uint32_t CRC32_ProcessBuffer(uint32_t ulCrc, const void *arg_pBuffer, unsigned int nBuffer);
