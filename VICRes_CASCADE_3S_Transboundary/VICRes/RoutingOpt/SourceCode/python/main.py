from sockets import Driver, InterfaceSocket, Status
import numpy as np
import time

n_reach=463
numres=45
sh = (1,n_reach+numres-1)
sha=(1,numres-1)
dtype = np.float32

def do_something(A):
    B = A[0,n_reach:n_reach+numres-1]
    return B*2

def main():
    server = InterfaceSocket()
    server.open()
    #time.sleep(240)
    client, address = server.server.accept()
    # client.settimeout(server.timeout)
    driver = Driver(client)
    print(" @SOCKET:   Client asked for connection from " + str(address) + ". Now hand-shaking.")

    A = np.zeros(sh, dtype)
    B = np.zeros(sha,dtype)
    # while(1):
    
    for i in range(2*732):
        stat = driver.get_status()
        if stat == Status.Up | Status.NeedsInit:
            driver.initialise()
            print('Driver initialised.')
        if stat == Status.Up | Status.HasData:
            #print('siamo qui e non va piu')
            A = driver.get_data(A)
            print('Data received: \n {}'.format(A))
            B = do_something(A)
        elif stat == Status.Up | Status.Ready:
            driver.send_data(B)
            print('Data sent: \n {}'.format(B))


if __name__ == '__main__':
    main()