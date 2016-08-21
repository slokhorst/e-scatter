#!/usr/bin/env python3
import struct

class record:
    def __init__(self, x,y,z,dx,dy,dz,K,xi,yi):
        self.x,self.y,self.z = x,y,z
        self.dx,self.dy,self.dz = dx,dy,dz
        self.K = K
        self.xi,self.yi = xi,yi

    def write(self, file):
        file.write(struct.pack('7f', self.x,self.y,self.z,self.dx,self.dy,self.dz,self.K))
        file.write(struct.pack('2i', self.xi,self.yi))
        file.flush()

    @classmethod
    def read(self, file):
        x, y, z, dx, dy, dz, K = struct.unpack('7f', file.read(7*4))
        xi, yi = struct.unpack('2i', file.read(2*4))
        return record(x,y,z,dx,dy,dz,K,xi,yi)

def read_all(filename):
    records = []
    with open(filename, "rb") as f:
        while f.read(1):
            f.seek(-1, 1)
            records.append(record.read(f))
    return records