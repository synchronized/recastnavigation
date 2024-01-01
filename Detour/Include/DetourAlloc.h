//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef DETOURALLOCATOR_H
#define DETOURALLOCATOR_H

#include <stddef.h>

/// Provides hint values to the memory allocator on how long the
/// memory is expected to be used.
/// 向内存分配器提供有关内存预计使用多长时间的提示值。
enum dtAllocHint
{
	///< Memory persist after a function call.
	///< 函数调用后内存仍然存在。
	DT_ALLOC_PERM,
	///< Memory used temporarily within a function.
	///< 函数内临时使用的内存。
	DT_ALLOC_TEMP
};

/// A memory allocation function.
///  @param[in]		size			The size, in bytes of memory, to allocate.
///  @param[in]		rcAllocHint	A hint to the allocator on how long the memory is expected to be in use.
///  @return A pointer to the beginning of the allocated memory block, or null if the allocation failed.
/// @see dtAllocSetCustom
/// 内存分配函数。
///  @param[in]		size          要分配的内存大小（以字节为单位）。
///  @param[in]		rcAllocHint	  向分配器提示内存预计使用多长时间。
///  @return 指向已分配内存块开头的指针，如果分配失败则为 null。
typedef void* (dtAllocFunc)(size_t size, dtAllocHint hint);

/// A memory deallocation function.
///  @param[in]		ptr		A pointer to a memory block previously allocated using #dtAllocFunc.
/// @see dtAllocSetCustom
/// 内存释放函数。
///  @param[in]		ptr		指向先前使用#dtAllocFunc分配的内存块的指针。
typedef void (dtFreeFunc)(void* ptr);

/// Sets the base custom allocation functions to be used by Detour.
///  @param[in]		allocFunc	The memory allocation function to be used by #dtAlloc
///  @param[in]		freeFunc	The memory de-allocation function to be used by #dtFree
/// 设置 Detour 使用的基本自定义分配函数。
///  @param[in]		allocFunc	#dtAlloc要使用的内存分配函数
///  @param[in]		freeFunc	#dtFree 使用的内存释放函数
void dtAllocSetCustom(dtAllocFunc *allocFunc, dtFreeFunc *freeFunc);

/// Allocates a memory block.
///  @param[in]		size	The size, in bytes of memory, to allocate.
///  @param[in]		hint	A hint to the allocator on how long the memory is expected to be in use.
///  @return A pointer to the beginning of the allocated memory block, or null if the allocation failed.
/// @see dtFree
/// 分配一个内存块。
///  @param[in]		size	要分配的内存大小（以字节为单位）。
///  @param[in]		hint	向分配器提示内存预计使用多长时间。
///  @return 指向已分配内存块开头的指针，如果分配失败则为 null。
void* dtAlloc(size_t size, dtAllocHint hint);

/// Deallocates a memory block.
///  @param[in]		ptr		A pointer to a memory block previously allocated using #dtAlloc.
/// @see dtAlloc
/// 释放内存块。
///  @param[in]		ptr		指向先前使用#dtAlloc分配的内存块的指针。
void dtFree(void* ptr);

#endif
