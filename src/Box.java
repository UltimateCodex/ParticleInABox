import java.awt.Color;
public class Box {

	public static void main(String[] args) {
		
	//Give Inputs
        int	   N  = 1024;		//Particle number
        double dt = 0.001;		//Time interval
        double T  = 100000;		//Run Time
		double G  = 10000;			//Self gravity
		double g  = 0;			//Downward gravity
		double k  = 0;		//Electric force
		double u  = 0;			//Magnetic force	c = sqrt(k/u)
		double Tw = 0;			//Wall Temperature
			
        //Set Background and Behavior
		StdDraw.setCanvasSize(500,500);
		boolean Wall = true;	
		boolean Well = false;	
		boolean Wrap = false;
		boolean Fuse = false;
		boolean Push = false;
			
		//Create Array of Particles	
		double A = 6;		//Configuration
		Particle[] P = new Particle[N];
//		Particle(x-position, y-position, x-velocity, y-velocity, mass, radius, elasticity, charge, image)
		for(int n = 0; n < N; n++){
			if(A == 0){P[n] = new Particle(Math.random(), Math.random(), 0, 0, 1, 0.005, 0, 1, "blueball.png");}					
			if(A == 1){P[n] = new Particle(Math.random(), Math.random(), (2*Math.random() - 1), (2*Math.random() - 1), 1, 0.01, 1, 1, "blueball.png");}
			if(A == 2){P[n] = new Particle((n / (int) Math.sqrt(N))/Math.sqrt(N) + 1/(2*Math.sqrt(N)), (n % Math.sqrt(N))/(Math.sqrt(N)) + 1/(2*Math.sqrt(N)), 0, 0, 1, 0.005, 1, 2*(n%2)-1, "blueball.png");}
			if(A == 3){P[n] = new Particle((n / (int) Math.sqrt(N))/Math.sqrt(N) + 1/(2*Math.sqrt(N)), (n % Math.sqrt(N))/(Math.sqrt(N)) + 1/(2*Math.sqrt(N)), 2*Math.random() - 1, 2*Math.random() - 1, 1, 0.005, 1, 1, "blueball.png");}
			if(A == 4){P[n] = new Particle((n / (int) Math.sqrt(N/2))/Math.sqrt(2*N) + 1/(2*Math.sqrt(2*N)), (2*n % Math.sqrt(2*N))/(2*Math.sqrt(2*N)), 0, 0, 0.0078125, 0.0078125, 1, 1, "blueball.png");}
			if(A == 5){P[n] = new Particle(0.5, 1.0*n/N, 0, 1, 1, 0.005, 0, 1, "blueball.png");}
			if(A == 6){P[n] = new Particle(Math.random(), Math.random(), 0.0001*(2*Math.random() - 1), 0.0001*(2*Math.random() - 1), Math.pow(10, 4*Math.random() - 10), -1, 0, 0, null);}
		}
		String form = "speed";
		
//		P[0] = new Particle(0.5, 0.5, 0, 0, 0, 0.05, -1, 10000, "blueball.png");
//		P[1] = new Particle(0.35, 0.5, 0, 0, 1, 0.05, 0, 1000, "blueball.png");
//		P[2] = new Particle(0.65, 0.5, 0, 0, -1, 0.05, 0, 1000, "blueball.png");

//		for(int n = 1; n < N; n++) {
//			P[n] = new Particle(0.2  + 0.01*(n%5), 1.5 - 0.01*(n/10), 0, -1, 0.00001, 0.005, -1, -0.00001, "blueball.png");
//		}
	
	//	P[0] = new Particle(0.199999, 0.53, 0.0, 0, 20, 0.01, 1, 0, "redballsmall.png");
	//	P[1] = new Particle(0.200001, 0.50, 0.0, 0, 40, 0.02, 1, 0, "greenballsmall.png");
	//	P[2] = new Particle(0.200000, 0.44, 0.0, 0, 80, 0.04, 1, 0, "blueball.png");
/*		
		P[3] = new Particle(0.399999, 0.53, 0.0, 0, 20, 0.01, 0.1, 0, "redballsmall.png");
		P[4] = new Particle(0.400001, 0.50, 0.0, 0, 40, 0.02, 1, 0, "greenballsmall.png");
		P[5] = new Particle(0.400000, 0.44, 0.0, 0, 80, 0.04, 1, 0, "blueball.png");
		
		P[6] = new Particle(0.599999, 0.53, 0.0, 0, 20, 0.01, 1, 0, "redballsmall.png");
		P[7] = new Particle(0.600001, 0.50, 0.0, 0, 40, 0.02, 0.1, 0, "greenballsmall.png");
		P[8] = new Particle(0.600000, 0.44, 0.0, 0, 80, 0.04, 1, 0, "blueball.png");
		
		 P[9] = new Particle(0.799999, 0.53, 0.0, 0, 20, 0.01, 1, 0, "redballsmall.png");
		P[10] = new Particle(0.800001, 0.50, 0.0, 0, 40, 0.02, 1, 0, "greenballsmall.png");
		P[11] = new Particle(0.800000, 0.44, 0.0, 0, 80, 0.04, 0.1, 0, "blueball.png");		
*/		
//		P[2] = new Particle(0.5, 0.53,     0, 0, 20, 0.01, 0.9, 0, "blueball.png");
//		P[1] = new Particle(0.5, 0.50,     0, 0, 40, 0.02, 0.9, 0, "blueball.png");
//		P[0] = new Particle(0.5, 0.44, 0.001, 0, 80, 0.04, 0.9, 0, "blueball.png");
		
//		P[2] = new Particle(0.499999, 0.23, 0.0, 0, 20, 0.01, 0.9, 0, "redballsmall.png");
//		P[1] = new Particle(0.500001, 0.20, 0.0, 0, 40, 0.02, 0.9, 0, "greenballsmall.png");
//		P[0] = new Particle(0.500000, 0.14, 0.0, 0, 80, 0.04, 0.9, 0, "blueball.png");
		
				
		//Run Simulation
		for(double t = 0; t < T; t += dt){
			
			//System.out.println(t);
			
			StdDraw.clear(StdDraw.GRAY);		
			
			Double[] Fx = new Double[N]; 
			Double[] Fy = new Double[N]; 
			
			for(int n = 0; n < N; n++){
		//		Fy[n] = -P[n].m*1;
		//		Fx[0] = P[n].NetForceX(P);
				Fx[n] = P[n].NetForceX(P,G,k);
				Fy[n] = P[n].NetForceY(P,G,k,g);
				Fx[0] = 0.0;
				Fy[0] = 0.0;
			}
			
			for(int n = 0; n < N; n++){
			for(int m = 0; m < N; m++){
				if( Push) {P[n].BallPushOff(P[m]);}
				if(!Fuse) {P[n].BallCollide(P[m]);}
				if( Fuse) {P[n].BallCombine(P[m]);}
			}
				P[n].WallCollide(Tw,Wall,Well,Wrap);
			}

			for(int n = 0; n < N; n++){
				P[n].Move(dt, Fx[n], Fy[n]);
			}
	/*		
			double Mx = 0;
			double My = 0;
			double M  = 0;
			double Kx = 0;
			double Ky = 0;
			double K  = 0;
			
			double MinV = 1e10;
			double AvgV = 0;
			double MaxV = 0;
			
			double Tp = 0;
			double Pr = 0;	
			
			for(Particle p: P){
				
				Mx += p.m*p.Vx;
				My += p.m*p.Vx;
				M  += p.m*p.V();
				Kx += 0.5*p.m*p.Vx*p.Vx;
				Ky += 0.5*p.m*p.Vy*p.Vy;
				K  += 0.5*p.m*p.V()*p.V();
				
				AvgV += p.V()/N;
				if(p.V() < MinV){MinV = p.V();}
				if(p.V() > MaxV){MaxV = p.V();}
				
				Tp += p.m*p.V()*p.V()/(3*N);				//k = 1
				Pr += p.WallCollide()/(4*t);
				
				
				final Color Temp = new Color(p.R(), p.G(), p.B());
		        StdDraw.setPenColor(Temp);
				StdDraw.filledCircle(p.Px, p.Py, p.r);
//				StdDraw.picture(p.Px, p.Py, p.img, 2*p.r, 2*p.r, p.o);
			}	
//			StdDraw.filledCircle(P[0].Px, P[0].Py, P[0].r);
			System.out.printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", 
							  Mx, My, M, Kx, Ky, K, MinV, AvgV, MaxV, Tp, Pr, N*Tp/Pr);
	*/
			
//	        StdDraw.text(0.8, 0.8, t);
			
			for(Particle p: P){
				if(form == "speed"){
					final Color Temp = new Color(p.R(), p.G(), p.B());
			        StdDraw.setPenColor(Temp);
			        StdDraw.filledCircle(p.Px, p.Py, p.r);
				}
				if(form == "image"){
					StdDraw.picture(p.Px, p.Py, p.img, 2*p.r, 2*p.r);
				}
				if(form == "spin"){
					StdDraw.picture(p.Px, p.Py, "colorwheel.png", 2*p.r, 2*p.r, p.o);
				}
			}
			for(int n = 0; n < N; n++){
				if(form == "split1"){
					if(n  < N/2){StdDraw.picture(P[n].Px, P[n].Py, "greenballsmall.png", 2*P[n].r, 2*P[n].r, P[n].o);}
					if(n >= N/2){StdDraw.picture(P[n].Px, P[n].Py, "blueballsmall.png", 2*P[n].r, 2*P[n].r, P[n].o);}
				}
				if(form == "split2"){
					if(n%2 == 0){StdDraw.picture(P[n].Px, P[n].Py, "greenballsmall.png", 2*P[n].r, 2*P[n].r, P[n].o);}
					if(n%2 == 1){StdDraw.picture(P[n].Px, P[n].Py, "blueballsmall.png", 2*P[n].r, 2*P[n].r, P[n].o);}
				}
			}
			StdDraw.show(0);
			
		}

	}
}
