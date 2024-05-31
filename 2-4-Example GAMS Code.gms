sets
         e index for Temperature                 /1*5/
         j index for Supplier                    /1*4/
         z index for Parts                       /1*2/
         l index for Light                       /1*5/;

scalars
         a  retailer's ordering cost             /46/
         b  vendor's ordering cost               /81/
         d  annual Demand rate                   /1420/
         f  product waste factor                 /0.4/
         h  retailer's holding cost              /8/
         m  vehicle capacity                     /114/
         r  vendor's holding cost                /59/
         s  warehouse carbon dioxide per product  /0.66/
         T  carbon tax coefficient                     /13/
         v  transportation carbon dioxide per product /0.15/
         Cw water consumption cost                     /6/
         wc water consumption in production           /20/
         ww water consumption in warehouse            /26/
         Clq     light cost per product                  /5/
         Ceq     temperature cost per product            /1/
         k     The cost of waste removal               /2664/
         ;

parameters
   Cl(l) light cost             /1 5,2 15,3 25,4 35,5 50/
   Ce(e) temperature cost       /1 40,2 70,3 100,4 130,5 160/

   AA(j) The jth supplier ordering cost         /1 100,2 101,3 108,4 104/
   tt(j) Transportation cost of jth supplier(vehicle)/1 99,2 83,3 84,4 91/
   w(z)  Occupied space by zth part(vehicle)     /1 2,2 4/;

table
         C(z,j)  purchasing cost of the z^{th} item from j^{th} supplier
         1          2          3          4
1        32         25         22         29
2        20         23         39         26;

table
         O(z,j)  supplier capacity
         1           2           3           4
1        686         597         650         469
2        403         560         228         348;

variables
         ZZ;

binary variables
         X(j)
         XX(z,j)

         U(l)   Light levels
         G(e)   Environment temperature levels ;

         U.l(l)=0;
         G.l(e)=0;

integer variables
         q      The retailer's order quantity
         n      Per cycle replenishment frequency of retailer
         P(z,j) Purchased order quantity of zth part from jth supplier;

         q.lo=1;
         q.up=1000;
         n.lo=1;
         P.up(z,j)=o(z,j);

Equations
         ObjectiveFunction
         coDemand(z)
         coCap(z,j)
         cobinaryX(j)
         cobinary(z,j)
         co1n
         co2U
         co3G   ;
         a=1.3*a;

ObjectiveFunction ..   ZZ =e= (h*q/2) + (b+(n*a))*(d/(n*q))

                            + ((r*q)*(n+1)/2)*(d/(n*q))

                            + ( T*((v*D) + s*( (q/2) + (q*(n+1)/2)*(D/(n*q)))) )

                            +( Cw*( ( (wc+ww)*(q*(n+1)/2))*(d/(n*q)) ) )

                            +  k*f*(1/(((Clq*(q*(n+1)/2)) + sum(l,U(l)*Cl(l)))

                            + ((Ceq*(q*(n+1)/2)) + sum(e,G(e)*Ce(e)))) ) *(d/(n*q))

                            + (      sum((z,j),C(z,j)*P(z,j))  + sum(j,AA(j)*X(j))

                            +   sum( j, tt(j) * ceil( sum(z, w(z)*P(z,j)) /m ) )     ) *(d/(n*q)) ;


         coDemand(z)       .. sum(j, P(z,j)) =e= n*q;
         coCap(z,j)        .. P(z,j) =l= o(z,j);
         cobinaryX(j)      .. X(j) =e= ifthen( sum(z,P(z,j))>0 ,1,0 );
         cobinary(z,j)     .. XX(z,j) =e= ifthen(P(z,j)>0 ,1,0 );
         co1n              .. n =l= 53;
         co2U              .. sum(l,U(l)) =e= 1 ;
         co3G              .. sum(e,G(e)) =e= 1 ;

model supply /all/;
option limrow=100, limcol=1000;

*option reslim=3600;
Option decimals=0;
solve supply using MINLP minimizing ZZ;
display ZZ.l, q.l, n.l, U.l, G.l, P.l, X.l, XX.l ;
display supply.resusd;