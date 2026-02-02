// Copyright Hitbox Games, LLC. All Rights Reserved.

#include "Modules/ModuleManager.h"

/**
 * Module implementation for CompressedSparseRow plugin.
 *
 * This is a header-only library module that provides CSR graph data structures.
 * No runtime initialization or shutdown logic is required.
 */
class FCompressedSparseRowModule : public IModuleInterface
{
};

IMPLEMENT_MODULE(FCompressedSparseRowModule, CompressedSparseRow)
