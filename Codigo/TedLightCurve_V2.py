#Ted Leandro @ UFRN june 1th 2017
#Lets make a automatic light curve analysis tool from scratch for Kepler's Stars

#we begin by importing the needed packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import kplr
import math
from astropy.stats import LombScargle
from numpy.polynomial import Chebyshev as T
import scipy.signal as signal
from gatspy.periodic import LombScargleFast

#definicoes de funcoes
#-##############################################################################
#-##############################################################################

def get_object(Object,ID,Qtin,Qtfin,method):

    client = kplr.API()
    if Object == "planet":
        box = client.planet(ID)

    if Object == "star":
        box = client.star(ID)

    lcs = box.get_light_curves()
    lcs = lcs[Qtin:Qtfin]

    # Loop over the datasets and read in the data.
    time, flux_a, flux_c, err_a, err_b = [], [], [], [], []
    s = 0
    for lc in lcs:
        with lc.open() as f:
            # coleto os dados de cada arquivo.
            hdu_data = f[1].data
            # aqui entra o tempo, fluxo e erro originais de cada quarter
            time = np.append(time, hdu_data["time"])
            flux_a = np.append(flux_a, hdu_data["pdcsap_flux"])
            err_a  = np.append(err_a,  hdu_data["pdcsap_flux_err"])
            if s == 0:
                cadencia = (max(time)-min(time))/len(time)
                s = 1

            # aqui entra a normalizacao e concatenacao para o caso de divisao
            flux_div, err_div = normalize(hdu_data["pdcsap_flux"], hdu_data["pdcsap_flux_err"], method)
            flux_c = np.append(flux_c, flux_div)
            err_b =  np.append(err_b, err_div)
    return time, flux_a, flux_c, err_a, err_b


def normalize(flux,err,method):
    """ Normaliza o fluxo de entrada pelo metodo de subtracao ou divisao
        se escolhido o metodo de divisao ("div"), ira normalizar em 1, e se escolhido
        o metodo de subtracao ("sub") ira normalizar em 1e6."""
    onde = 1e8# 671688 #onde quero normalizar no caso de subtracao (soma)

    if method == "div":
        mediaF = np.nanmedian(flux)
        mediaE = np.nanmedian(err)
        flux = flux/mediaF
        err = err/mediaE - 1

    if method == "sub":
        mediaF = np.nanmedian(flux)
        flux = flux + (onde - mediaF)

    return flux, err

def reduction (time_nan, flux_nan, err, cadencia, factor=5, norm=1):

    #agora vou encontrar as gaps entre os quarters e outras gaps menores e preencher com um determinado ruido
    faltante_time = [] #esta lista guardara o tempo que falta nos gaps
    i = 1       #comecarei a partir do 1 para pegar o segundo comeco (que eh do segundo quarter) pois o primeiro nao interessa
    # neste for, irei conferir cada intervalo entre os pontos fotometricos para conferir se ha alguma pausa maior do que um determinado tempo
    #de observacao, se houver, entao devo preencher essa gap com um determinado ruido
    for i in range (len(time_nan)-1):
        teste = time_nan[i]-time_nan[i-1] #este teste confere o valor de tempo entre dois pontos fotometricos
        # se o intervalo entre os pontos for maior do que a cadencia definida da observacao, entao entro no if para preencher esse espaco
        if teste > 0.0405:
            valor = time_nan[i]-time_nan[i-1] #guarda o intervalo de tempo que precisarei preencher
            falta_time = np.linspace(time_nan[i-1],time_nan[i],(valor/cadencia)) #cria uma lista que contem os horarios que deveriam haver observacao baseado na cadencia predefinida
            faltante_time = np.append(faltante_time, falta_time) #inclui esta lista na minha lista de intervalos das gaps

    faltante_flux = np.zeros(len(faltante_time))+norm # crio uma lista com a mesma quantidade de pontos que o tempo que falto nos gaps normalizados em norm (padrao = 1)
    faltante_err = np.zeros(len(faltante_time))
    faltante_flux = np.random.normal(faltante_flux,np.std(flux_nan)/factor) #aplico alguma randomicidade nesse fluxo
    time_nan = np.append(time_nan, faltante_time) #incluo o tempo sintetico dos gaps no tempo original
    flux_nan = np.append(flux_nan, faltante_flux) #incluo o fluxo sintetico dos gaps no fluxo original
    err_nan = np.append(err, faltante_err)

    time_nan, flux_nan = zip(*sorted(zip(time_nan, flux_nan))) #organizo de forma temporal minha lista total de tempo e fluxo
    time_nan = np.array(time_nan)   #garanto o aspecto de array para o tempo
    flux_nan = np.array(flux_nan) #garanto o aspecto de array para o fluxo flux

    return time_nan, flux_nan, err_nan


def detrending(time,flux, err, pl):

    fluxTrue = np.isfinite(flux) #confere se em um index tem um numero ou nao. retorna True ou False
    index = np.where(fluxTrue == True)[0] #indexs que possuem numeros (True)
    flux_nan = flux[index] #todos os valores de fluxo sem os "nan"s
    time = time[index]          # o tempo que corresponde ao fluxo sem "nan"
    err_nan = err[index]
    time_nan = (time - min(time))   # normaliza o inicio do tempo em zero dias
    # aqui eu coloco em ordem numerica temporal a curva de luz total
    time_nan, flux_nan = zip(*sorted(zip(time_nan, flux_nan)))
    time_nan = np.array(time_nan)   #garanto o aspecto de array para o tempo
    flux_nan = np.array(flux_nan) #garanto o aspecto de array para o fluxo flux


    p = T.fit(time_nan, flux_nan, pl) # faco aqui uma fitagem da curva de forma polinomial com o polinomio definido
    flux_model = p(time_nan)        #aplico a minha fitagem com o tempo de observacao
    flux_detrended = (flux_nan-flux_model)   # faco a subtracao do real com o modelo
    flux_detrended = flux_detrended + 1       # como isso vai pra zero, normalizo de volta em 1

    return time_nan, flux_nan, flux_model, flux_detrended, err_nan

def detrending2(time,flux, err, pl):

    p = T.fit(time, flux, pl) # faco aqui uma fitagem da curva de forma polinomial com o polinomio definido
    flux_model = p(time)        #aplico a minha fitagem com o tempo de observacao
    flux_detrended = (flux-flux_model)   # faco a subtracao do real com o modelo
    flux_detrended = flux_detrended + 1       # como isso vai pra zero, normalizo de volta em 1

    return time, flux, flux_model, flux_detrended, err


def binagem (time_done, flux_done, err, period, nbins):
    foldtime = time_done/period
    foldtime = foldtime % 1
    width = 1.0/nbins #tamanho de cada bin

    bins = np.zeros(nbins)
    weights = np.zeros(nbins)

    for i in range(len(flux_done)):
        n = int(foldtime[i] / width)
        weight = err[i]**-2.0
        bins[n] += flux_done[i]*weight
        weights[n] += weight

    bins /= weights

    binErr = np.sqrt(1.0/(weights))
    binEdges = (np.arange(nbins)*width)
    binEdges = np.linspace(0,period*24,nbins)

    binEdges = binEdges
    plt.plot(binEdges, bins,"og")
    #plt.errorbar(binEdges,bins,yerr=binErr,linestyle='none',marker='o')  # plot binned lightcurve
    plt.show()

#-##############################################################################
#-##############################################################################
Object = "planet"
ID = "kepler-666b"
Qtin = 3
Qtfin = 18
method = "div"
detrendy = 1000
preencher = 50
binn = 4000

# aqui posso chamar as funcoes e fazer o que eu quero

time, flux_original, flux, err_original, err = get_object(Object, ID, Qtin, Qtfin, method)
time_nan, flux_nan, flux_model, flux_detrended, err = detrending(time,flux, err, detrendy)
time_done, flux_done, err = reduction(time_nan,flux_detrended,err,preencher,factor=4)

#time_done, flux_done, flux_model2, flux_detrended2, err = detrending2(time_done,flux_done, err, 300)
#flux_done = flux_detrended2

#ploto essas reducoes
plt.subplot(3, 1, 1)
plt.title("curva de luz original, normalizada e detrended")
plt.plot(time, flux_original, 'kx')
plt.subplot(3, 1, 2)
plt.plot(time_nan,flux_nan,'kx')
plt.plot(time_nan,flux_model,'b-',lw=1.5)
plt.ylim(0.99,1.01)
plt.subplot(3, 1, 3)
plt.plot(time_nan, flux_detrended, 'kx')
plt.ylim(0.99,1.01)
plt.show()

# essa frequencia de busca mostra todo o espectro de potencia.
"""frequencia de Nyquist para gigantes = 24.4512 para inverso de dia"""
frequency = np.linspace(0.01, 1, 1000)
power = LombScargle(time_done, flux_done).power(frequency)
#frequency, power = LombScargle(time_done, flux_done).autopower()
plt.plot((1/frequency),(power))
plt.xlim(0,50)
#plt.xscale("log")
plt.show()


"""trabalando aqui em fitar o melhor range para o Prot """
#
"""
#ja descobri que as frequencias de busca ficam entre 0.004567 e 0.015222
frequency = np.linspace(0.004567, 0.015222, 100000)
power = LombScargle(time_done, flux_done).power(frequency)
#frequency, power = LombScargle(time_done, flux_done).autopower()
plt.loglog(1/frequency, power)
#plt.xscale("log")
plt.show()
"""

#-#########################################

period = 16.09 #periodo do kepler-131b do guilherme
#period = 25.5169
#period = 63.14
#period = 94.9   #periodo da binaria 5006817 do bruno
#period = 143.191 #Prot da estrela 5307747 do bruno
#period = 133.12   #Prot da estrela 5006817 do bruno
#period = 3.52187 #periodo do kepker-8b para quartes > 10 e para o comeco eh 3.5261
#period = 4.8851
#period = 3.54839
#period = 7.9725 # periodo do kepler-210b
#period = 2.4532 # periodo do kepler-210c
#period = 4.4988301  #periodo do planeta 1 do kepler-666
period = 2.25061   #periodo do segundo planeta que EU descobri no kepler-666
#period = 34.59
#period = 229
#period = 3.5225
#period = 20.8851
#period = 16.87

foldtime = time_done/period
foldtime = foldtime % 1

plt.subplot(2, 1, 1)
plt.ylim(0.99,1.01)
plt.plot(time_done, flux_done,'rx')
plt.subplot(2, 1, 2)
plt.plot(foldtime, flux_done, 'kx')
plt.ylim(0.99,1.01)
plt.show()

binagem(time_done,flux_done,err,period,binn)

#-#################################










