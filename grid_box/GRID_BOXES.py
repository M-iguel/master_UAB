import sys
box_size=int(sys.argv[1])
grid_size=int(sys.argv[2])

def box(box_size,grid_size):
    print(("+" + "-" * (box_size-2))*grid_size + "+")#imprimo la primera línea de +--+--+...
    for t in range(grid_size):#quiero hacer una caja tantas veces como sea el grid_size
        i=0
        j=box_size-4 #j será los espacios entre \ y / y empieza en box_size-4 porque ya he hecho 2 lados y he añadido \ y / (el resto serán espacios)
        for v in range(box_size//2-1):#haré la primera mitad de las filas menos la primera, que ya la he hecho
            print (grid_size*("|"+i*" "+"\\"+j*" "+"/"+i*" "),"|",sep="")#imprimo las líneas en las que \ va antes que / añadiendo espacios progresivamente crecientes a la izquierda y decrecientes en el medio de ambas. Al final añado el lado del final | sin separación
            i+=1#la i marca los espacios antes de \ y va creciendo
            j-=2#la j marca los espacios ente \ y / y va decreciendo
        if box_size%2 != 0:#si el resto box_size/2 es distinto de 0 es q el nº es impar
            for p in range(box_size%2):#en ese caso solo quiero una iteración (hago un range q vaya del 0 al 1, que es el resto)
                print (grid_size*("|"+i*" "+"X"+j*" "+i*" "),"|",sep="")#imprimo la X en el medio
            k=1 #empiezo en k=1 pq acabo de añadir una línea en el medio de la caja para poner la X
            l=box_size-4
            for v in range(box_size//2-1):#esto es para hacer la segunda mitad menos el final, que lo haré después
                print (grid_size*("|"+(l-i)*" "+"/"+k*" "+"\\"+(l-i)*" "),"|",sep="")#ahora añado primero / y luego \ y los espacios los empiezo a hacer donde los dejó la otra (por eso resto i)
                k+=2#ahora los espacios del medio van creciendo
                l-=1#los espacios de los lados van decreciendo (estamos en la mitad inferior del cuadrado)
        else: #si el número es par
            k=0 #empiezo en k=0 porque no he añadido nunguna línea
            l=box_size-4
            for v in range(box_size//2-1):#esto es para hacer la segunda mitad menos el final, que lo haré después
                print (grid_size*("|"+(l-(i-1))*" "+"/"+k*" "+"\\"+(l-(i-1))*" "),"|",sep="")#ahora añado primero / y luego \ y los espacios los empiezo a hacer donde los dejó la otra
                k+=2
                l-=1
        #última línea es igual que la primera
        print("+",end='')#imprimo un + sin que cambie de línea
        for n in range(box_size-2):#quiero tantas líneas como mi ancho-2 (que son las esquinas)
            print ("-",end='')#voy imprimiendo - en la misma línea
        print(("+" + "-" * (box_size-2)) * (grid_size-1) + "+")#imprimo la última línea

box(box_size,grid_size)

