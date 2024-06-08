#durée est le temps en secondes durant lequel des mesures seront effectuées
#nbdonnées est le nombre de données à enregistrer

def Mesures(durée):
    import time
    t0 = time.time()                           #stocke le nombre de secondes écoulées depuis le 01/01/1970
    import serial                              #importation du module serial pour arduino
    serialArduino = serial.Serial('COM5',9600) #identification du port arduino
    Ltof=[]                                    #liste des temps en microsecondes
    Lang=[]                                    #liste des angles en degrés
    a=int(print("Entrer le nombre de valeurs à enregistrer"))
    while type(a)!=type(1) or a>2 or a<1:
        a=int(print("Entrer le nombre de valeurs à enregistrer"))
    serialArduino.write("démarrage")           #instruction envoyée à arduino initiant son démarrage
    while t<=durée:
        données=serial_port.readline().split()
        try:
            tof=données[0]
            ang=données[1]
            Ltof.append(tof)
            Lang.append(ang)
        except:
            pass
        t=time.time()-t0
    serialArduino.write("arrêt")           #instruction envoyée à arduino entrainant son arrêt
    serial_port.close()
    if a==1:
        return(Ltof)
    else:
        return
