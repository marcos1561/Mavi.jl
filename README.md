# Mavi

[![Build Status](https://github.com/marcos1561/Mavi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcos1561/Mavi.jl/actions/workflows/CI.yml?query=branch%3Amain)

Mavi é um motor de dinâmica de partículas (_Particle Dynamic Engine_).

# Interface visual
O Mavi possui uma interface visual (feita inteiramente com o [Makie](https://docs.makie.org/v0.21/)) cujo objetivo é servir de ferramenta de depuração visual para o sistema sendo explorado. A estrutura da UI possui essencialmente dois elementos:

1. Gráfico onde as partículas do sistema são renderizadas em tempo real.
2. Painel de informações sobre o estado do sistema e execução do programa.

<!-- ![alt text](docs/images/ui_components.png "Title") -->
<img src="docs/images/ui_components.png" alt="Componentes da UI" width="500"/>

## Animando o sistema
Dado que um função `step!(system, int_cfg)` já está construída para algum `system`, podemos animar o processo de integração da seguinte forma

```julia
using Mavi.Visualization

# Criando o sistema
system = ...

animate(system, step!)
```

É possível configurar aspectos da animação passando uma instância de `AnimationCfg` em `animate`. O seguinte exemplo anima o sistema com o fps setado para 30, executando 15 passos temporais a cada frame da animação 

```julia
using Mavi.Visualization

# Criando o sistema
system = ...

anim_cfg = AnimationCfg(
    fps = 30,
    num_steps_per_frame = 15,
)

animate(system, step!, anim_cfg)
```

## Estendendo o painel de informações
É possível injetar informações customizáveis no painel de informações. Fazemos isso setando o campo `custom_items` de `DefaultInfoUICfg`, que por sua vez é um campo de `AnimationCfg`. `custom_items` é uma função que deve retornar as informações adicionais que serão mostradas no painel de informações. Para informações mais detalhadas de sua assinatura, consulte sua documentação em [info_ui.jl](src/gui/info_ui.jl).  
O seguinte exemplo utiliza um sistema já definido no Mavi e adiciona a informação da posição da primeira partícula no painel de informações

```julia
using Mavi
using Mavi.Configs
using Mavi.Visualization
using Printf

system = System(
    state = State{Float64}(
        pos = [[1 2 3]; [1 2 3]],
        vel = [[1 1 0]; [-1 0 2]],
    ),
    space_cfg = RectangleCfg(length=4, height=4),
    dynamic_cfg = LenJonesCfg(sigma=1, epsilon=0.1),
    int_cfg = IntCfg(dt=0.01),
)

function get_pos(system, _)
    pos = system.state.pos[:, 1]
    pos_formatted = @sprintf("(%.3f, %.3f)", pos[1], pos[2])
    return [("pos_1", pos_formatted)]
end

anim_cfg = AnimationCfg(
    info_cfg = DefaultInfoUICfg(
        custom_items = get_pos
    )
)

animate(system, Mavi.Integration.step!, anim_cfg)
```
