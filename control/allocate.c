/**
 *
 * @file allocate.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void *plasma_shared_alloc(plasma_context_t *plasma, size_t size, int type)
{
    void *memptr;

    size *= plasma_element_size(type);
    if (size <= 0)
        return NULL;
  //if (posix_memalign(&memptr, STANDARD_PAGE_SIZE, size) != 0) {
    if ((memptr = malloc(size)) == NULL) {
        plasma_error("plasma_shared_alloc", "posix_memalign() failed");
        return NULL;
    }
	if ( plasma->runtime == PLASMA_OMPSS) {
		#pragma omp register([size]memptr)
//		printf("shared_alloc::memptr: %p[%d]\n", memptr, size);
	}
    return memptr;
}

/***************************************************************************//**
 *
 **/
void plasma_shared_free(plasma_context_t *plasma, void *ptr)
{
    if (ptr == NULL)    // somewhat redundant - free() does the same
        return;
	if ( plasma->runtime != PLASMA_OMPSS) {
		free(ptr);
    }
}

/***************************************************************************//**
 *
 **/
void *plasma_private_alloc(plasma_context_t *plasma, size_t size, int type)
{
    void *memptr;

    size *= plasma_element_size(type);
    if (size <= 0)
        return NULL;
  //if (posix_memalign(&memptr, CACHE_LINE_SIZE, size) != 0) {
    if ((memptr = malloc(size)) == NULL) {
        plasma_error("plasma_private_alloc", "posix_memalign() failed");
        return NULL;
    }
	if ( plasma->runtime == PLASMA_OMPSS) {
		#pragma omp register([size]memptr)
//		printf("private_alloc::memptr: %p[%d]\n", memptr, size);
	}
    return memptr;
}

/***************************************************************************//**
 *
 **/
void plasma_private_free(plasma_context_t *plasma, void *ptr)
{
    if (ptr == NULL)    // somewhat redundant - free() does the same
        return;
	if ( plasma->runtime != PLASMA_OMPSS) {
		free(ptr);
    }
}
