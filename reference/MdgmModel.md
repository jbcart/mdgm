# MDGM Model

MDGM Model

MDGM Model

## Details

R6 class wrapping a C++ MDGM model. Use
[`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md)
to create instances.

## Methods

### Public methods

- [`MdgmModel$new()`](#method-MdgmModel-new)

- [`MdgmModel$has_emission()`](#method-MdgmModel-has_emission)

- [`MdgmModel$nvertices()`](#method-MdgmModel-nvertices)

- [`MdgmModel$ncolors()`](#method-MdgmModel-ncolors)

- [`MdgmModel$emission_type()`](#method-MdgmModel-emission_type)

- [`MdgmModel$get_ptr()`](#method-MdgmModel-get_ptr)

------------------------------------------------------------------------

### Method `new()`

Create a new MdgmModel. Use
[`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md)
instead.

#### Usage

    MdgmModel$new(ptr)

#### Arguments

- `ptr`:

  External pointer to a C++ Model object.

------------------------------------------------------------------------

### Method `has_emission()`

Check if the model has an emission layer.

#### Usage

    MdgmModel$has_emission()

#### Returns

Logical.

------------------------------------------------------------------------

### Method `nvertices()`

Get the number of vertices.

#### Usage

    MdgmModel$nvertices()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `ncolors()`

Get the number of colors (categories).

#### Usage

    MdgmModel$ncolors()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `emission_type()`

Get the emission family type.

#### Usage

    MdgmModel$emission_type()

#### Returns

Character string (`"bernoulli"`, `"gaussian"`, `"poisson"`) or `NULL`
for standalone models.

------------------------------------------------------------------------

### Method `get_ptr()`

Get the internal C++ pointer. For internal use only.

#### Usage

    MdgmModel$get_ptr()

#### Returns

External pointer.
