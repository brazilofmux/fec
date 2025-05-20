/* This code is public domain -- Will Hartung 4/9/09 */
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

size_t getline(char** lineptr, size_t* n, FILE* stream)
{
    if (lineptr == nullptr || stream == nullptr || n == nullptr)
    {
        return -1;
    }
    char *bufptr = *lineptr;
    size_t size = *n;

    int c = fgetc(stream);
    if (c == EOF)
    {
        return -1;
    }
    if (bufptr == nullptr)
    {
        bufptr = new char[128];
        if (bufptr == nullptr)
        {
            return -1;
        }
        size = 128;
    }
    size_t i = 0;
    while (c != EOF)
    {
        if (size <= i+1)
        {
	        const size_t resize = size + 128;
            auto bufptr2 = new char[resize];
            if (bufptr2 == nullptr)
            {
                delete[] bufptr;
                bufptr = nullptr;
                return -1;
            }
            memcpy(bufptr2, bufptr, size);
            delete[] bufptr;
            bufptr = bufptr2;
            bufptr2 = nullptr;
            size = resize;
        }
        if (c == '\n' || c == '\r')
        {
            break;
        }
        bufptr[i++] = static_cast<char>(c);
        c = fgetc(stream);
    }

    bufptr[i++] = '\0';
    *lineptr = bufptr;
    *n = size;

    return i;
}