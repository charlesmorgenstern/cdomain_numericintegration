using Plots
using Roots



function getpaths(n1,n2,n3,n4) #generate four paths that define the c domain
#paths defined in counterclockwise orientation for Green's Thm

paths=Dict();

t=collect(range(pi/4,2*pi,n1)); #outer circle of radius 3
x=3*cos.(t)
y=3*sin.(t)
p1=Matrix{Float64}(undef,n1,2)
p1[:,1]=x
p1[:,2]=y
paths["p1"]=p1;


x=collect(range(3,1,n2));  #horizontal line
p2=Matrix{Float64}(undef,n2,2)
p2[:,1]=x
p2[:,2].=0.0
paths["p2"]=p2;

t=t[end:-1:1] #inner circle of radius one
t=collect(range(2*pi,pi/4,n3));
x=cos.(t)
y=sin.(t)
p3=Matrix{Float64}(undef,n3,2)
p3[:,1]=x
p3[:,2]=y
paths["p3"]=p3;


x=collect(range(cos(pi/4),3*cos(pi/4),n4)); #line of slope 1
p4=Matrix{Float64}(undef,n4,2)
p4[:,1]=x
p4[:,2]=x
paths["p4"]=p4;

return paths

end

################################################################
################################################################


function plotdomain() #plots the domain legend shows indexing of paths

paths=getpaths(100,100,100,100)
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]
fig=plot(p1[:,1],p1[:,2])
plot!(p2[:,1],p2[:,2])
plot!(p3[:,1],p3[:,2])
plot!(p4[:,1],p4[:,2])
return fig

end

################################################################
################################################################


function greensintegral(n1,n2,n3,n4)  #ni points for the ith path
#numeric approx of integral for area using Green's Thm
#both double integral and boundary can be varified analytically by hand 
#to be 7pi. We make no assumptions of what the curves are.
#We only assume the points are equally spaced in some parameter t.
#We let t range from 0 to 1 for each curve.

paths=getpaths(n1,n2,n3,n4) #load paths with ni points 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line 
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

int1=0.0 #line integral for outer circle
for i=2:n1
x=p1[i,1]
y=p1[i,2]
dx=(p1[i,1]-p1[i-1,1])/h1 #first order finite difference approx of derivatives
dy=(p1[i,2]-p1[i-1,2])/h1
int1+=x*dy-y*dx
end 
int1*=h1

int2=0.0 #line integral for horizontal line
for i=2:n2
x=p2[i,1]
y=p2[i,2]
dx=(p2[i,1]-p2[i-1,1])/h2 #first order finite difference approx of derivatives
dy=(p2[i,2]-p2[i-1,2])/h2
int2+=x*dy-y*dx
end 
int2*=h2

int3=0.0 #line integral for inner circle
for i=2:n3
x=p3[i,1]
y=p3[i,2]
dx=(p3[i,1]-p3[i-1,1])/h3 #first order finite difference approx of derivatives
dy=(p3[i,2]-p3[i-1,2])/h3
int3+=x*dy-y*dx
end 
int3*=h3


int4=0.0 #line integral for line w/ slope 1
for i=2:n4
x=p4[i,1]
y=p4[i,2]
dx=(p4[i,1]-p4[i-1,1])/h4 #first order finite difference approx of derivatives
dy=(p4[i,2]-p4[i-1,2])/h4
int4+=x*dy-y*dx
end 
int4*=h4

int=.5*(int1+int2+int3+int4)


exactarea=7*pi

display("Exact area of c domain is:")
display(exactarea)
display("Numeric approx of Green's Thm integral of c domain is:")
display(int)
display("Relative error is:")
err=abs((exactarea-int)/exactarea)
display(err)
end

###############################################################
###############################################################


function montecarlo_cdom(n) #monte carlo approx of area integral
#plots random points on domain labelling which points are inside domain

ptin=(-2,0) #generate random points in [-4,4]x[-4,4]
pts=(rand(Float64,(n,2)))
pts=8*(pts.-.5)


n1=100 #define number of points for the paths
n2=100
n3=100
n4=100

paths=getpaths(n1,n2,n3,n4) #generate paths
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]


fig=plot(p1[:,1],p1[:,2]) #initialize plot with domain
plot!(p2[:,1],p2[:,2])
plot!(p3[:,1],p3[:,2])
plot!(p4[:,1],p4[:,2])



ptsin=0 #start counting how many points are in domain

for j=1:n #loop through random points

#create line from random point to point inside domain
g(x)=((pts[j,2]-ptin[2])/(pts[j,1]-ptin[1]))*(x-ptin[1])+ptin[2]
k=0 #start counting intersections with boundary of domain

#check for intersections with outer circle
for i=2:n1 #connect points defining path with simple lines
f(x)=((p1[i,2]-p1[i-1,2])/(p1[i,1]-p1[i-1,1]))*(x-p1[i-1,1])+p1[i-1,2]
h(x)=f(x)-g(x) #get difference of two lines
v=find_zeros(h,p1[i-1],p1[i]) #roots of difference are intersections
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1])) 
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n2 #check intersections with horizontal line
f(x)=((p2[i,2]-p2[i-1,2])/(p2[i,1]-p2[i-1,1]))*(x-p2[i-1,1])+p2[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p2[i-1],p2[i])
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n3 #check intersections with inner circle
f(x)=((p3[i,2]-p3[i-1,2])/(p3[i,1]-p3[i-1,1]))*(x-p3[i-1,1])+p3[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p3[i-1],p3[i])
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n4 #check intersections with line of slope 1
f(x)=((p4[i,2]-p4[i-1,2])/(p4[i,1]-p4[i-1,1]))*(x-p4[i-1,1])+p4[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p4[i-1],p4[i])
cnt=length(v)
if cnt>=1&&((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end

if iseven(k) #if number of intersections is even point is in domain
ptsin+=1 #count point
scatter!([pts[j,1]],[pts[j,2]], mc=:red,legend = false)
else
scatter!([pts[j,1]],[pts[j,2]], mc=:blue,legend = false)
end
#display(k)
end
ratio=ptsin/n
area=64*ratio
display("approximate area of domain is")
display(area)
display("exact area of domain is")
display(7*pi)
display("relative error is")
err=abs(area-7*pi)/(7*pi)
display(err)
return fig
end

###############################################################
###############################################################


function montecarlo_noplot(n) #same as montecarlo_cdom but without plotting

ptin=(-2,0)
pts=(rand(Float64,(n,2)))
pts=8*(pts.-.5)


n1=100
n2=100
n3=100
n4=100

paths=getpaths(n1,n2,n3,n4)
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]





ptsin=0

for j=1:n

g(x)=((pts[j,2]-ptin[2])/(pts[j,1]-ptin[1]))*(x-ptin[1])+ptin[2]
k=0

for i=2:n1
f(x)=((p1[i,2]-p1[i-1,2])/(p1[i,1]-p1[i-1,1]))*(x-p1[i-1,1])+p1[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p1[i-1],p1[i])
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1])) 
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n2
f(x)=((p2[i,2]-p2[i-1,2])/(p2[i,1]-p2[i-1,1]))*(x-p2[i-1,1])+p2[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p2[i-1],p2[i])
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n3
f(x)=((p3[i,2]-p3[i-1,2])/(p3[i,1]-p3[i-1,1]))*(x-p3[i-1,1])+p3[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p3[i-1],p3[i])
cnt=length(v)
if cnt>=1 && ((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end


for i=2:n4
f(x)=((p4[i,2]-p4[i-1,2])/(p4[i,1]-p4[i-1,1]))*(x-p4[i-1,1])+p4[i-1,2]
h(x)=f(x)-g(x)
v=find_zeros(h,p4[i-1],p4[i])
cnt=length(v)
if cnt>=1&&((pts[j,1]<ptin[1] && v[1]<ptin[1] && v[1]>pts[j,1]) ||(pts[j,1]>ptin[1] && v[1]>ptin[1] && v[1]<pts[j,1]))
if cnt>1
display("Error: possible tangent")
display(pts[j,:])
end
k+=cnt
end
end

if iseven(k)
ptsin+=1
end
end
ratio=ptsin/n
area=64*ratio
display(area)
display("exact area of domain is")
display(7*pi)
display("relative error is")
err=abs(area-7*pi)/(7*pi)
display(err)
end



