// Copyright Hitbox Games, LLC. All Rights Reserved.

#include "Modules/ModuleManager.h"

/**
 * Module implementation for CompactSparseRow plugin.
 *
 * This is a header-only library module that provides CSR graph data structures.
 * No runtime initialization or shutdown logic is required.
 */
class FCompactSparseRowModule : public IModuleInterface
{
};

IMPLEMENT_MODULE(FCompactSparseRowModule, CompactSparseRow)
