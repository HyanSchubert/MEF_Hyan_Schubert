using LinearAlgebra

include("material.jl")
include("quad_enriquecido.jl")
include("apoios.jl")
include("global.jl")

function main()

    # Entrada de dados na mão!

    nnos = 8
    ne = 3

    coord = [0.0   0.0   ;
             1.0   0.0   ;
             2.0   0.0   ;
             3.0   0.0   ;
             0.0   5.0   ;
             1.0   5.0   ;
             2.0   5.0   ;
             3.0   5.0   ]                  


    conectividades = [1 2 6 5;
                      2 3 7 6;
                      3 4 8 7]

    VE =     ones(ne)
    Vnuxy = (1/3)*ones(ne)

    # Espessura
    esp = 1.0

    # Apoios (cond. de contorno essenciais)
    #        no gl valor
    apoios = [1 1 0.0;
              1 2 0.0;
              2 2 0.0;
              3 2 0.0;
              4 2 0.0;
              5 1 0.0]

    # Forças (cond. de contorno naturais)
    #         no gl valor
    P = 1.0
    forcas = [4  1  P;
              8  1  P]

    # Monta a matriz global do problema
    K = Global(ne,nnos,conectividades,coord, VE, Vnuxy, esp)

    # Monta o vetor de forças globais concentradas
    F = Forca_global(nnos,forcas)

    # Modifica o sistema pela aplicação das condições de contorno
    # homogêneas
    # Aplica_CC_homo!(apoios,K,F)
    KA, FA = Aplica_CC_lagrange(nnos,apoios,K,F)

    # Soluciona o sistema de equações, obtendo os deslocamentos
    # KU = F
    UA = KA\FA

    # Só os deslocamentos que nos interessam
    U = UA[1:2*nnos]

    # Retorna U
    return U

end
