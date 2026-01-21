# Rezbits

Rezbits is a compression, decompression, etc. library that currently
supports only DEFLATE decompression.

The library's name comes from [this](https://youtu.be/qk_Dwspn6jc)
Metroid Prime 2 enemy.

## Usage

This is an `#include` only library. The idea is for you to copy the
source folder, rename it to Rezbits, and add it to your project. Then
you would include the header files that you want like this:
`#include "Rezbits/deflate.hpp"`