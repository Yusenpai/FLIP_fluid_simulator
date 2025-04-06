# FLIP fluid simulator

A FLIP fluid simulator ported to C++, based on the implementation from [Ten Minute Physics](https://www.youtube.com/watch?v=XmzBREkK8kY&t=506s).

## Install

The demo use SFML graphic library. More infomation at https://www.sfml-dev.org. Install the SFML version 2.6.2. Build:

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```