// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
/** A result as JSON.parse gives it. */
export type Json = null | boolean | number | string | Json[] | { [k: string]: Json };
/** Runs the function `name` (its R name) on its named arguments; throws an Error with the reason when it refuses them. */
export function jrun(name: string, args?: Record<string, unknown> | string): Json;
/** The functions and their argument names ("?" marks an optional one). */
export function jfns(): Record<string, string[]>;
/** Every function as a method, by its R name. */
declare const crimlab: Readonly<Record<string, (args?: Record<string, unknown> | string) => Json>>;
export default crimlab;
