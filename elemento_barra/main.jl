using LinearAlgebra

include("pre.jl")
include("barra.jl")
include("global.jl")
include("apoios.jl")

# Função main que calcula o problema de barra por CC homogeneas
function main()

    # Entrada de dados na mão!
    ne = 3
    nnos = 3
    conectividades = [1 2;
                      2 3;
                      1 3]

    coord = [0.0 0.0;
             1.0 0.0;
             1.0 1.0]                  

    VE = 210E9*ones(ne)
    VA = 1E-4*ones(ne)

    # Pré-processamento
    VL,Vtheta = Pre_processa(ne,coord,conectividades)

    # Condições de contorno essenciais
    #        nó gll valor
    apoios = [1 1 0;
              1 2 0;
              2 2 0]
    
    # Condições de contorno naturais
    #        nó gll valor
    forcas = [2 1 1000.0;
              3 1 1000.0;
              3 2 1000.0]

    # Monta a matriz global do problema
    K = Rigidez_global(ne,nnos,conectividades, VE, VA, VL, Vtheta)

    # Monta vetor de força global do problema
    F = Forca_global(nnos, forcas)
    
    # Aplica_CC_lagrange
    KA, FA = Aplica_CC_lagrange(nnos,apoios,K,F)

    # Soluciona o sistema de equações, obtendo os deslocamentos
    # KU = F
    UA = KA\FA
   
    # Só os deslocamentos que nos interessam
    UL = UA[1:2*nnos]

    # Aplica Rearranjo
    Rearranjo!(nnos,apoios,K,F)

    # Soluciona o sistema de equações, obtendo os deslocamentos
    # KU = F
    U = K\F

    return UL

end #main


