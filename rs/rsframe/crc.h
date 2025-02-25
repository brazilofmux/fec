// Copyright (C) 1993-1996 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

#ifndef CRC_H
#define CRC_H

UINT32 CRC32_ProcessBuffer
(
    UINT32         ulCrc,
    const void    *arg_pBuffer,
    unsigned int   nBuffer
);

#endif // CRC_H
