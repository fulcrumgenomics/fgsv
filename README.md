# fgsv

Tools to find evidence for structural variation.

## Building & Testing

This repo uses [mill](https://com-lihaoyi.github.io/mill/mill/Intro_to_Mill.html) as it's build system.

To run unit tests:

```console
./mill _.test
```

To build an executable JAR at `./jars/fgsv.jar`:

```console
./mill _.deployLocal
```
