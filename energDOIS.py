# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy
import math


'''
SEGUNDA PARAMETRIZACAO

Grafico De Energia de Ligacao por Numero de Nucleons de
acordo com o Modelo Bethet-Wiezsache

O programa busca implementar a funcao do modelo para
energia de ligacao nuclear um termo de cada vez, e 
gerar um grafico com um termo a mais por vez para
permitir a visualizacao do efeito dos diversos termos.



'''


#Valores dos parametros
a_v = 15.5
a_s = 16.8
a_c = 0.72
a_i = 23.0



def delta(N,Z):
    '''
    Essa funcao implementa o termo delta da funcao, que
    tem o valor de acordo com a paridade dos numeros de
    nucleons
    '''
    A = N+Z #Numero total de Nucleons
    delt = 0 #inicializacao do valor de delta
    
    if (N+Z)%2 == 0: #Se o caso for par-par ou impar-impar 
        if N%2 == 0: #Se for par-par
           delt = 34*A**(-3/4.)

           #No caso par-par eh adicionado um termo

        else: #Se nao eh impar-impar
            delt = -34*A**(-3/4.)

            #No caso impar impar ocorre um decremento

    #Se o if não for acionado é o caso impar-par,
    #e delta = 0

    
    return delt


def energ_vol(N,Z):
    #Retorna o termo de volume
    
    return a_v*(N+Z)

def energ_sup(N,Z):
    #Retorna o Termo anterior + superficie
    A = N+Z
    return a_v*A - a_s*A**(2/3.)


def energ_col(N,Z):
    #Retorna o termo anterior + termo coulombiano
    A = N+Z
    div = A**(1/3.)
    add = a_c*Z*(Z-1)/div
    return energ_sup(N,Z) -add

def energ_sim(N,Z):
    #Retorna Termo anterior + termo de simetria
    A = N+Z
    return energ_col(N,Z) - a_i*(N-Z)*(N-Z)/A

def energ_del(N,Z):
    #Retorna termo anterior + termo delta
    return energ_sim(N,Z) + delta(N,Z)



def read_column(filename):
    '''
    Essa função cria uma lista com uma coluna de inteiros
    num arquivo de texto
    '''
    
    f = open(filename,"r") #Abre arquivo
    text = f.read().replace("\r","").replace(" ","")
    #LÊ e remove espaços e termos desnecessários
    
    split_text = text.split("\n")
    #Divide a lista de acordo com quebras de linha

    column = []
    for elem in split_text:
        column.append(int(elem))

    #Transforma cada elemento da lista em inteiros e retorna ela
    return column


def create_nucleons(proton_list, neutron_list):
    '''
    Dado dois arquivo de texto com uma coluna de número de
    nucleons ele cria uma duas listas, uma com o número de
    protons e outras com número de neutrons

    No caso vou usar os arquivos 'protons.dat' e 'neutrons.dat',
    com os valores dos isótopos mais comunsx


    '''
    
    filename1 = proton_list
    filename2 = neutron_list
    #Nomes dos arquivos
    Protons = read_column(filename1)
    Neutrons = read_column(filename2)
    #Chama a funçao anterior neles
    
    if len(Protons)==len(Neutrons):
        print('As listas tem o mesmo tamanho')
        return Protons,Neutrons

    else:
        print('As listas NAAOO tem o mesmo tamanho')
        return Protons, Neutrons
    #Verifica se as listas tem o mesmo tamanho
    


def plt_energ(Proton_list, Neutron_list):
    '''
    Dado arquivos com número de nucleons, cria listas para
    o número de protons e neutrons e calcula para cada caso
    da função de energia por número barionico.

    Em seguida gera um gráfico disso
    '''
    Protons, Neutrons = create_nucleons(Proton_list,Neutron_list)    
    #Cria as listas com a função acima

    en_vol = []
    en_sup = []
    en_col = []
    en_sim = []
    en_del = []
    A=[]
    #Cria listas vazias para guardas o valores
    for elem in Protons:
        z = elem
        n = Neutrons[elem-1]
        a = n+z
        A.append(n+z)
        en_vol.append(energ_vol(n,z)/a)
        en_sup.append(energ_sup(n,z)/a)
        en_col.append(energ_col(n,z)/a)
        en_sim.append(energ_sim(n,z)/a)
        en_del.append(energ_del(n,z)/a)

        #Calcula a energia/A e salva na lista respectiva


    f = plt.figure() #inicia uma figura

    area = 8 #area dos pontos
    alp = 0.5 #Transparencia dos pontos
    
    plt.scatter(A, en_vol,s=area,alpha =alp, label = 'Energia Vol.')
    plt.scatter(A, en_sup,s=area,alpha =alp, label = '-Superficial')
    plt.scatter(A, en_col,s=area/2,alpha =alp, label = '-Coloumbiano')
    plt.scatter(A, en_sim,s=area/2,alpha =alp, label = '-Simetrico')
    plt.scatter(A,en_del,s=area/2,alpha =alp, label = '+Resto')
    #Printando os 5 casos
    
    plt.ylim(ymin = 3, ymax = 18) #Limite do eixo y
    
    plt.xlabel('A')
    plt.ylabel('Energia/A')
    plt.title('Isotopo Mais Comum: Caso 2')

    plt.legend(loc="lower center")
    
    f.savefig('graf-caso2.pdf')


plt_energ('protons.dat','neutrons.dat')

