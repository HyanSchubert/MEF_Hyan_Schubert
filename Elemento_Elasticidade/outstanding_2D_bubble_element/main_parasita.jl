using LinearAlgebra

include("material.jl")
include("quad_enriquecido.jl")
include("apoios.jl")
include("global.jl")

function exemplo_cis_parasita()

    # Entrada de dados na mão!

    nnos = 22
    ne = 10
    L=2.5
    coord = [0.0   0.0   ; #1
             0.25  0.0   ; #2
             0.5   0.0   ; #3
             0.75  0.0   ; #4
             1.0   0.0   ; #5
             1.25  0.0   ; #6
             1.5   0.0   ; #7
             1.75  0.0   ; #8
             2.0   0.0   ; #9
             2.25  0.0   ; #10
             2.5   0.0   ; #11
             2.5   0.2   ; #12
             2.25  0.2   ; #13
             2.0   0.2   ; #14
             1.75  0.2   ; #15  
             1.5   0.2   ; #16
             1.25  0.2   ; #17
             1.0   0.2   ; #18
             0.75  0.2   ; #19
             0.5   0.2   ; #20
             0.25  0.2   ; #21
             0.0   0.2   ] #22                


    conectividades = [1 2 21 22;
                      2 3 20 21;
                      3 4 19 20;
                      4 5 18 19;
                      5 6 17 18;
                      6 7 16 17;
                      7 8 15 16;
                      8 9 14 15;
                      9 10 13 14;
                      10 11 12 13]  

    E=210E9
    VE = E*ones(ne)
    Vnuxy = (1/3)*ones(ne)

    # Espessura
    h = 0.05

    # Apoios (cond. de contorno essenciais)
    #        no gl valor
    apoios = [1 1 0.0;
              1 2 0.0;
              22 1 0.0;
              22 2 0.0]

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

    # Retorna U
    return Un_max, Ua_max

end
