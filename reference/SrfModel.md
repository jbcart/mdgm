# Spatial Random Field Model

Spatial Random Field Model

Spatial Random Field Model

## Details

R6 class wrapping a C++ spatial random field model. Use
[`srf_model()`](https://jbcart.github.io/mdgm/reference/srf_model.md) to
create instances.

## Methods

### Public methods

- [`SrfModel$new()`](#method-SrfModel-new)

- [`SrfModel$has_emission()`](#method-SrfModel-has_emission)

- [`SrfModel$nvertices()`](#method-SrfModel-nvertices)

- [`SrfModel$ncolors()`](#method-SrfModel-ncolors)

- [`SrfModel$emission_type()`](#method-SrfModel-emission_type)

- [`SrfModel$model_type()`](#method-SrfModel-model_type)

- [`SrfModel$get_ptr()`](#method-SrfModel-get_ptr)

------------------------------------------------------------------------

### Method `new()`

Create a new SrfModel. Use
[`srf_model()`](https://jbcart.github.io/mdgm/reference/srf_model.md)
instead.

#### Usage

    SrfModel$new(ptr, model_type = "mdgm")

#### Arguments

- `ptr`:

  External pointer to a C++ Model object.

- `model_type`:

  Character string: `"mdgm"` or `"mrf"`.

------------------------------------------------------------------------

### Method `has_emission()`

Check if the model has an emission layer.

#### Usage

    SrfModel$has_emission()

#### Returns

Logical.

------------------------------------------------------------------------

### Method `nvertices()`

Get the number of vertices.

#### Usage

    SrfModel$nvertices()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `ncolors()`

Get the number of colors (categories).

#### Usage

    SrfModel$ncolors()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `emission_type()`

Get the emission family type.

#### Usage

    SrfModel$emission_type()

#### Returns

Character string (`"bernoulli"`, `"gaussian"`, `"poisson"`) or `NULL`
for standalone models.

------------------------------------------------------------------------

### Method `model_type()`

Get the spatial model type.

#### Usage

    SrfModel$model_type()

#### Returns

Character string: `"mdgm"` or `"mrf"`.

------------------------------------------------------------------------

### Method `get_ptr()`

Get the internal C++ pointer. For internal use only.

#### Usage

    SrfModel$get_ptr()

#### Returns

External pointer.

## Methods

### Public methods

- [`MdgmModel$new()`](#method-SrfModel-new)

- [`MdgmModel$has_emission()`](#method-SrfModel-has_emission)

- [`MdgmModel$nvertices()`](#method-SrfModel-nvertices)

- [`MdgmModel$ncolors()`](#method-SrfModel-ncolors)

- [`MdgmModel$emission_type()`](#method-SrfModel-emission_type)

- [`MdgmModel$model_type()`](#method-SrfModel-model_type)

- [`MdgmModel$get_ptr()`](#method-SrfModel-get_ptr)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    MdgmModel$new(ptr, model_type = "mdgm")

------------------------------------------------------------------------

### Method `has_emission()`

#### Usage

    MdgmModel$has_emission()

------------------------------------------------------------------------

### Method `nvertices()`

#### Usage

    MdgmModel$nvertices()

------------------------------------------------------------------------

### Method `ncolors()`

#### Usage

    MdgmModel$ncolors()

------------------------------------------------------------------------

### Method `emission_type()`

#### Usage

    MdgmModel$emission_type()

------------------------------------------------------------------------

### Method `model_type()`

#### Usage

    MdgmModel$model_type()

------------------------------------------------------------------------

### Method `get_ptr()`

#### Usage

    MdgmModel$get_ptr()
