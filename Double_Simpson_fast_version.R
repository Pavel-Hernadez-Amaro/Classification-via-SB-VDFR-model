
# NA_ind=NULL
# for (ind_i in 1:dim(Real_X[[j]])[1]) {
#   for (ind_j in 1:dim(Real_X[[j]])[2]) {
#     if(is.na(Real_X[[j]][ind_i,ind_j])){
#       NA_ind=sort(c(NA_ind,(dim(Real_X[[j]])[1]*(ind_j-1))+ind_i))
#     }}}
# 

n=2*sub # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
m=2*sub

simp_w=rep(1,m+1)
even=seq(2,m+1-1,2)
odd=seq(3,m+1-1,2)

simp_w[even]=2
simp_w[odd]=4

Sim_w_x=(h/3)*diag(simp_w)
Sim_w_y=(HX/3)*diag(simp_w)

h=(x_b-x_a)/n
HX=(y_b-y_a)/m # if y_b and y_a are functions of x
# this line should go inside the next for loop (the int_i lopp or outer loop)

x = x_a + (0:n)*h
y = y_a + (0:m)*HX

na_x=which(x>=n_1[1]) # ESTO ESTA PARA UNA SOLA PERSONA HAY QUE AGREGAR LAS COORDENADAS FALTANTES DE CADA SUPERFICIES (MATRICES O LISTAS)
na_y=which(y<n_2[1])

L=matrix(h/3,nrow = c1*c2, ncol = c1*c2)

for (index in 1:N) {
  
  fx2=matrix(0,nrow = m+1, ncol = c2)
  
  fx2[1:length(na_y),]=spline.des(knots_y, y[na_y], 3+1, 0*y[na_y])$design
  
  
for (int_i in 1:length(x)) {

  if (int_i%in%na_x) {
    next
  }

    print(c("int_i =",int_i))


    #### first point
    fx1=spline.des(knots_x, x[int_i], 3+1, 0*x[int_i])$design
    
    fx=(kronecker(fx2,fx1))
    
    Fx=t(fx) %*% Sim_w_y %*% fx
 
    if (int_i==1 || int_i==length(x)) {
      
      L= L + Fx
      
    }else{
    if (int_i%%2==0) {
      L = L + 2*Fx
    }else{
      L = L + 4*Fx
    }}
    

  } # for in int_i
  

} # for in index

