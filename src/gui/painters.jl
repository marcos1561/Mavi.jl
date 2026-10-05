abstract type PainterCfg end

abstract type Painter end
abstract type AbstractPalettePainter <: Painter end

get_painter(cfg::PainterCfg, system) = get_painter(cfg, system.type, system)
get_painter(cfg::PainterCfg, sys_type, system) = get_painter(cfg, system.state, system.dynamic_cfg, system)
function get_painter(cfg::PainterCfg, state, dynamic_cfg, system)
    throw(ArgumentError("Not implemented for $(typeof(cfg))"))
end

function update_painter!(painter::Painter, system) end

# ===
# Palette Painter
# ===

@kwdef struct PaletteTypes{E} <: AbstractVector{Int}
    num_types::Int
    types::Vector{Int}
    extra::E=nothing
end

with_extra(p::PaletteTypes, extra) = PaletteTypes(
    num_types=p.num_types,
    types=p.types,
    extra=extra
)

Base.size(types::PaletteTypes) = size(types.types)
Base.IndexStyle(::Type{PaletteTypes}) = IndexLinear()
Base.getindex(types::PaletteTypes, i::Int) = types.types[i]
function Base.setindex!(types::PaletteTypes, value, i::Int)
    types.types[i] = value
    return types
end

abstract type TypesMode end

always_update_types(types_mode) = true
get_types(types_mode, state, dynamic_cfg, system) = get_types(types_mode, state)
update_types!(types_mode::TypesMode, painter::Painter, system) = update_types!(types_mode, painter.types, system)
update_types!(types_mode::TypesMode, types::PaletteTypes, system) = false

function get_types(types_mode::Vector, state) 
    # num_t, num_p = length(types_mode), length(state.pos)
    # if num_t != num_p
    #     error("length(types)=$(num_t), but it should be equal to the maximum number of particles $(num_p).")
    # end

    PaletteTypes(num_types=maximum(types_mode), types=types_mode)
end

function get_types(types_mode::Int, state; fill_value=1) 
    PaletteTypes(num_types=types_mode, types=fill(fill_value, length(state.pos)))
end

@kwdef struct EntityTypes <: TypesMode 
    gap::Int = 0
end
always_update_types(types_mode::EntityTypes) = false

function get_types(types_mode::EntityTypes, state)
    types = collect(1:length(state.pos)) .+ types_mode.gap
    PaletteTypes(num_types=types[end], types=types)
end

function update_types!(types_mode::EntityTypes, types::PaletteTypes, system)
    types .= vec(1:length(system.state.pos)) .+ types_mode.gap
    return true
end

@kwdef struct RandomTypes{R} <: TypesMode
    num_types::Int
    rng::R=Random.default_rng()
end
function get_types(types_mode::RandomTypes, state) 
    PaletteTypes(num_types=types_mode.num_types, types=rand(types_mode.rng, 1:types_mode.num_types, length(state.pos)))
end

abstract type TypesUpdater end

function update_types!(updater::TypesUpdater, painter, system) 
    throw(ArgumentError("Not implemented for $(typeof(updater))"))
end

struct FuncTypesUpdater{F} <: TypesUpdater 
    func::F
end
update_types!(updater::FuncTypesUpdater, painter, system) = updater.func(painter, system)


struct DynamicTypes{U<:TypesUpdater, T, E} <: TypesMode 
    updater::U
    types::T
    extra::E
end
function DynamicTypes(;updater, types, extra=nothing)
    if types isa Symbol
        types = symbol_to_types(types)
    end

    if updater isa Function
        updater = FuncTypesUpdater(updater)
    end

    DynamicTypes(updater, types, extra)
end
get_types(types_mode::DynamicTypes, state) = get_types(types_mode.types, state)
update_types!(types_mode::DynamicTypes, painter::Painter, system) = update_types!(types_mode.updater, painter, system)


abstract type ScalarUpdater end
struct NoUpdate <: ScalarUpdater end

update_scalar!(updater::NoUpdate, buffer, system) = false
update_scalar!(updater::Function, buffer, system) = updater(buffer, system)

@kwdef struct FakeContinuosTypes{N, U} <: TypesMode
    num_colors::Int=256
    gap::Int=0
    normalizer::N=identity
    updater::U=NoUpdate()
end

function get_types(types_mode::FakeContinuosTypes, state)
    types = get_types(types_mode.num_colors + types_mode.gap, state)
    types = with_extra(types, fill(0.0, length(types)))
end
function update_types!(types_mode::FakeContinuosTypes, types::PaletteTypes, system)
    buffer = types.extra
    to_update = update_scalar!(types_mode.updater, buffer, system)
    to_update = to_update === nothing ? true : to_update
    if to_update
        for (idx, s) in enumerate(buffer)
            x = types_mode.normalizer(s)
            types[idx] = round(Int, (types_mode.num_colors - 1) * x + 1) + types_mode.gap
        end
    end
end


struct ManyTypes{T} <: TypesMode 
    list::T
end
ManyTypes(; kwargs...) = ManyTypes((; kwargs...))

Base.getindex(types::ManyTypes, i::Union{Int, Symbol}) = types.list[i]

function get_types(types_mode::ManyTypes, state) 
    types_list = map(t -> get_types(t, state), types_mode.list)

    PaletteTypes(
        num_types=maximum(t -> t.num_types, values(types_list)),
        types=copy(first(types_list).types),
        extra=types_list,
    )
end

function symbol_to_types(symbol)
    types = nothing
    if symbol == :entities
        types = EntityTypes()
    else
        error("Types `$(type)` is not a valid type.")
    end

    return types
end


abstract type PaletteMode end

get_palette(palette_mode, state, dynamic_cfg, types, system) = get_palette(palette_mode, types)

function get_palette(palette_mode::Vector, types) 
    palette = []
    for c in palette_mode
        @show c
        push!(palette, RGBf(GLMakie.to_color(c)))
    end
    return palette
end

function get_palette(palette_mode::PaletteMode, types)
    throw(ArgumentError("Not implemented for $(typeof(palette_mode))"))
end

struct ManyPalettes <: PaletteMode
    list::Vector
end

get_palette(palettes::ManyPalettes, types) = vcat([get_palette(p, types) for p in palettes.list]...)

@kwdef struct RandomPalette{R} <: PaletteMode 
    length_offset::Int=0
    rng::R=Random.default_rng()
end

function get_palette(palette_mode::RandomPalette, types)
    length = mod1(types.num_types + palette_mode.length_offset, types.num_types)
    [rand(palette_mode.rng, RGBf) for _ in 1:length]
end

struct CmapPalette{C, S, R} <: PaletteMode 
    cmap::C
    sampler::S
    range::Tuple{Float64, Float64}
    length_offset::Int
    rng::R
end
function CmapPalette(; cmap, sampler=:random, range=(0.0, 1.0), length_offset=0, rng=Random.default_rng())
    a, b = Float64.(range)
    0.0 <= a <= b <= 1.0 ||
        throw(ArgumentError("range must satisfy 0 ≤ a ≤ b ≤ 1"))
    
    cmap = cmap isa Symbol ? colorschemes[cmap] : cmap

    CmapPalette(cmap, sampler, (a, b), length_offset, rng)
end

function get_palette(p::CmapPalette, types)
    sampler = p.sampler
    n = mod1(types.num_types + p.length_offset, types.num_types)
    a, b = p.range
    
    if sampler isa Symbol
        if sampler == :random
            sampler = (idx, rng) -> rand(rng)
        elseif sampler == :linear
            buffer = LinRange(0, 1, n)
            sampler = (idx, rng) -> buffer[idx]
        else
            throw(ArgumentError("sampler must be a function or :random or :linear; got $sampler"))
        end
    end
    
    map(1:n) do idx
        u = sampler(idx, p.rng)
        (u isa Real && 0 <= u <= 1) ||
            throw(ArgumentError(
                "sampler(rng) must return a finite real number in [0, 1]; got $u"
            ))

        get(p.cmap, a + u * (b-a))
    end
end

struct PalettePainterCfg{P, T} <: PainterCfg 
    palette::P
    types::T
end
function PalettePainterCfg(; palette=:random, types=:entities, rng=nothing)
    if rng === nothing
        rng = Random.default_rng()
    end

    if palette isa Symbol
        if palette === :random
            palette = RandomPalette(rng=rng)
        else
            cmap = nothing
            try
                cmap = colorschemes[palette]
            catch 
                error("Palette `$(palette)` is not a valid palette name, neither a valid colorsheme name.")
            end
            palette = CmapPalette(cmap=cmap, rng=rng)
        end 
    end

    if palette isa Vector && (palette[1] isa PaletteMode || palette[1] isa Vector)
        palette = ManyPalettes(palette)
    end

    if types isa Symbol
        types = symbol_to_types(types)
    end

    PalettePainterCfg(palette, types)
end
PalettePainterCfg(p; rng) = PalettePainterCfg(palette=p, rng=rng)
PalettePainterCfg(p::Nothing; rng) = PalettePainterCfg(rng=rng)
PalettePainterCfg(p::PalettePainterCfg; rng) = p


struct PalettePainter{P, T, C, E} <: AbstractPalettePainter 
    cfg::PalettePainterCfg{P, T}
    types::PaletteTypes{E}
    palette::Vector{C}
    colors::Vector{C}
    is_variable_number::Bool
end

function ScalarPainter(;updater, cmap=:viridis, normalizer=identity, num_colors=256)
    function func(painter, system)
        sg.update_types!(painter.cfg.types.types, painter, system)
        if paint_border
            border_color_type(painter, system)
        end
    end

    PalettePainterCfg(
        types=FakeContinuosTypes(
            num_colors=num_colors,
            updater=updater,
            normalizer=normalizer,
        ),
        palette=CmapPalette(cmap=cmap, sampler=:linear),
    )
end


function get_painter(cfg::PalettePainterCfg, state, dynamic_cfg, system)
    types = get_types(cfg.types, state, dynamic_cfg, system)
    palette = get_palette(cfg.palette, state, dynamic_cfg, types, system)
    colors = Vector{eltype(palette)}(undef, length(types.types))
    
    if types.num_types > length(palette)
        error("There are $(types.num_types) different types, but length of palette is $(length(palette)).")
    end

    # @show types.types[1:3]
    # @show palette[1:3]

    p = PalettePainter(cfg, types, palette, colors, is_variable_number(state))
    update_color_buffer!(p, get_particles_ids(system))
    return p
end

function update_color_buffer!(painter::PalettePainter, particle_ids)
    colors = painter.colors
    for (i, p_id) in enumerate(particle_ids)
        colors[i] = painter.palette[painter.types[p_id]]
    end
end

function update_painter!(painter::PalettePainter, system) 
    types_changed = false
    if always_update_types(painter.cfg.types)
        types_changed = update_types!(painter.cfg.types, painter, system)
        types_changed = types_changed === nothing ? true : types_changed
    end

    if types_changed || painter.is_variable_number
        update_color_buffer!(painter, get_particles_ids(system))
    end
end
