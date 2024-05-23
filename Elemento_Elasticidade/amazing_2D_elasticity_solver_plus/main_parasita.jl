using LinearAlgebra

include("material.jl")
include("quad.jl")
include("apoios.jl")
include("global.jl")

function main_parasita()

    # Entrada de dados na mão!

    nnos = 12
    ne = 5
    L=5
    coord = [0.0   0.0   ; #1
             1.0   0.0   ; #2
             2.0   0.0   ; #3
             3.0   0.0   ; #4
             4.0   0.0   ; #5
             5.0   0.0   ; #6
             5.0   1.0   ; #7
             4.0   1.0   ; #8
             3.0   1.0   ; #9
             2.0   1.0   ; #10
             1.0   1.0   ; #11
             0.0   1.0   ] #12                


    conectividades = [1 2 11 12;
                      2 3 10 11;
                      3 4 9 10;
                      4 5 8 9;
                      5 6 7 8;]

    E=210E9
    VE = E*ones(ne)
    Vnuxy = (1/3)*ones(ne)

    # Espessura
    h = 0.05

    # Apoios (cond. de contorno essenciais)
    #        no gl valor
    apoios = [1 1 0.0;
              1 2 0.0;
              12 1 0.0;
              12 2 0.0]

    # Forças (cond. de contorno naturais)
    #         no gl valor
    F = 1000
    forcas = [7  2  F]

    # Deslocamento máximo ANALÍTICO
    Ua_max = (F*L^3)/(3*E*(h^3/12))

    # Monta a matriz global do problema
    K = Global(ne,nnos,conectividades,coord, VE, Vnuxy, h)

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

    # Deflexão máxima NUMÉRICA
    Un_max=maximum(U)

    # Deslocamento máximo analítico
    #F=1000; L=5; E=210E9; h=0.05;
    #Utmax = (F*L^3)/(3*E*(h^3/12))

    # Retorna U
    return Un_max, Ua_max

end
