
public class Particle {
	
	public double Px;	//x position
	public double Py;	//y position
	public double Vx;	//x velocity
	public double Vy;	//y velocity
	public double m;	//mass
	public double r;	//radius
	public double e;	//elasticity
	public double q; 	//charge
	public double o = 0;	//orientation
	public double w = 0;	//angular velocity
	public double a = 1; 	//moment of Inertia
	public double u = 1;	//friction coefficient
	public String img;	//image
	
	public Particle(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, String p9){
		Px  = p1;
		Py  = p2;
		Vx  = p3;
		Vy  = p4;
		m   = p5;
		r   = p6;
		e   = p7;
		q   = p8;
		img = p9;
		
		if(p6 == -1) {r = Math.pow(m, 1.0/3);}
	}
	
	public Particle (Particle p){
		Px  = p.Px;
		Py  = p.Py;
		Vx  = p.Vx;
		Vy  = p.Vy;
		m   = p.m;
		r   = p.r;
		e   = p.e;
		q   = p.q;
		img = p.img;
	}
		
	public double V(){
		return Math.sqrt(Vx*Vx + Vy*Vy);
	}
	
	public double Dx(Particle p){
		return p.Px - Px;
	}
	
	public double Dy(Particle p){
		return p.Py - Py;
	}
	
	public double D(Particle p){
		return Math.sqrt(Dx(p)*Dx(p) + Dy(p)*Dy(p));
	}
	
	public double ForceN(Particle p, double G, double k){
		if(D(p) > p.r + r && (G != 0 || k != 0)){
			 return (G*p.m*m - k*p.q*q)/(D(p)*D(p));
		}
		else{return 0;}
	}
	
	public double ForceM(Particle p, double u){
		if(D(p) > p.r + r && u != 0){
			return (u*q*p.q*V()*p.V())/(D(p)*D(p));
		}
		else{return 0;}
	} 
	
	public double ForceX(Particle p, double G, double k){
		return ForceN(p,G,k)*Dx(p)/D(p);
	}
		
	public double ForceY(Particle p, double G, double k){
		return ForceN(p,G,k)*Dy(p)/D(p);
	}		
		
	public double NetForceX(Particle[] P, double G, double k){	
		double Fx = 0;
		for (Particle p: P){
			if(p.Px != Px || p.Py != Py){
				Fx += ForceX(p,G,k);
			}
		}	
		return Fx;
	}
	
	public double NetForceY(Particle[] P, double G, double k, double g){	
		double Fy = -m*g;
		for (Particle p: P){
			if(p.Px != Px || p.Py != Py){
				Fy += ForceY(p,G,k);
			}
		}
		return Fy;
	}
	
	public void Move(double dt, double Fx, double Fy){
		if(m != 0){
			Vx += dt*Fx/m;
			Vy += dt*Fy/m;
		}
			Px += dt*Vx;
			Py += dt*Vy;
			o  += dt*w;
	}		
	
	double P = 0;
	public double WallCollide(double T, boolean Wall, boolean Well, boolean Wrap){
		if(Wall){
			if(Px < r && Vx < 0 || Px > 1 - r && Vx > 0){
				Vx = -Vx + T;
				P += Math.abs(2*m*Vx);
			}
			if(Py < r && Vy < 0 || Py > 1 - r && Vy > 0){
				Vy = -Vy + T;
				P += Math.abs(2*m*Vy);
			}
		}
		if(Well){
			if(Px < r && Vx < 0 || Px > 1 - r && Vx > 0){
				Vx = -Vx + T;
			}
			if(Py < r && Vy < 0){
				Vy = -Vy + T;
			}
		}
		if(Wrap){
			if(Px < -r){Px = 1 + r;}
			if(Px > 1 + r){Px = -r;}
			if(Py < -r){Py = 1 + r;}
			if(Py > 1 + r){Py = -r;}
		}
		return P;
	}
	
	public void BallPushOff(Particle p){
		if(D(p) < p.r + r && m != 0){
			p.Px = Px + (p.r + r)*Dx(p)/D(p);
			p.Py = Py + (p.r + r)*Dy(p)/D(p); 
		}
	}
	public void BallCollide(Particle p){
		if(D(p) < p.r + r && Dx(p)*(p.Vx - Vx) + Dy(p)*(p.Vy - Vy) < 0){
						
			double E;
			if(p.e < e){E = p.e;} else{E = e;}
			
			double Vn1 = (  Vx*Dx(p) +   Vy*Dy(p))/D(p);
			double Vt1 = (  Vx*Dy(p) -   Vy*Dx(p))/D(p);
			double Vn2 = (p.Vx*Dx(p) + p.Vy*Dy(p))/D(p);			
			double Vt2 = (p.Vx*Dy(p) - p.Vy*Dx(p))/D(p);
			
			double Vfn1 = (m*Vn1 + p.m*Vn2 + p.m*E*(Vn2 - Vn1))/(m + p.m);
			double Vfn2 = (m*Vn1 + p.m*Vn2 +   m*E*(Vn1 - Vn2))/(m + p.m);
		
			double Vft1 = Vt1 + u*(1 + e)*(p.m/(m + p.m))*(Vn1 - Vn2);				
			double Vft2 = Vt2 + p.u*(1 + p.e)*(m/(m + p.m))*(Vn2 - Vn1);	
			
		//	  w = 360*(w - u*(1 + e)*(p.m/(m + p.m))*(Vn1 - Vn2)/(a*r));	
		//	p.w = 360*(p.w - p.u*(1 + p.e)*(m/(m + p.m))*(Vn2 - Vn1)/(p.a*p.r));	
			
	//		double Vft1 = ((1 - a*u)*Vt1 + a*(1 + u)*r*w)/(1 + a);
	//		double Vft2 = ((1 - p.a*p.u)*Vt2 + p.a*(1 + p.u)*p.r*p.w)/(1 + p.a);
			
	//		  w = 360*((1 + u)*Vt1 +   (a - u)*r*w)        /(r*(1 + a));
	//		p.w = 360*((1 + p.u)*Vt2 + (p.a - p.u)*p.r*p.w)/(p.r*(1 + p.a));
			
			  Vx = (Vfn1*Dx(p) + Vft1*Dy(p))/D(p);
			  Vy = (Vfn1*Dy(p) - Vft1*Dx(p))/D(p);
			p.Vx = (Vfn2*Dx(p) + Vft2*Dy(p))/D(p);
			p.Vy = (Vfn2*Dy(p) - Vft2*Dx(p))/D(p);
		}
	}
	
	public void BallCombine(Particle p){
		if(D(p) < p.r + r && Dx(p)*(p.Vx - Vx) + Dy(p)*(p.Vy - Vy) < 0 && m != 0){
			
			  Px = (m*Px + p.m*p.Px)/(m + p.m);
			  Py = (m*Py + p.m*p.Py)/(m + p.m);
			  Vx = (m*Vx + p.m*p.Vx)/(m + p.m);
			  Vy = (m*Vy + p.m*p.Vy)/(m + p.m);
		
			p.Px = 0.0;
			p.Py = 0.0;
			p.Vx = 0.0;
			p.Vy = 0.0;
								  
			  r  = r*Math.pow((m + p.m)/m, 1.0/3); 			
			p.r  = 0.0;	  
			
			  m += p.m;
			p.m  = 0.0;
		}
	}
	
	public int R(){int R = (int)(255*V());			if(R > 255){R = 255;}	if(R < 0){R = 0;}	return R;}
	public int G(){int G = (int)(255*(V() - 1));	if(G > 255){G = 255;}	if(G < 0){G = 0;}	return G;}
	public int B(){int B = (int)(255*(V() - 2));	if(B > 255){B = 255;}	if(B < 0){B = 0;}	return B;}
	
	//	R = (int)(255*p.V()/20);
	//	B = (int)(255*(1 - p.V()/20));	
		
}
		
