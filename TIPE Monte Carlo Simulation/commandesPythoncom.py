# -*- coding: utf-8 -*-
import serial

# class Arduino():
def __init__(self,port):
    self.ser = serial.Serial(port,baudrate=57600)
    c_recu = self.ser.read(1)
    while ord(c_recu)!=0:
        c_recu = self.ser.read(1)
    c_recu = self.ser.read(1)
    while ord(c_recu)!=255:
        c_recu = self.ser.read(1)
    c_recu = self.ser.read(1)
    while ord(c_recu)!=0:
        c_recu = self.ser.read(1)
    self.SEND_BUFFER = 100
    self.READ_SEND_BUFFER = 101
    self.N = 256 # taille des tableaux
        
def close(self):
    self.ser.close()
