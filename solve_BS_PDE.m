function V = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, right_cond, left_bound, right_bound, T, Mx, Mt)
	 
	 ## usage: V = solve_BS_PDE (term_cond, left_cond, right_cond,
  ## left_bound, right_bound, T, Mx, Mt)
	 ## 
	 ## 
	 
	 xmin = min(left_cond);
	 xmax= max(right_cond);
	 V = term_cond';
	 x = linspace(xmin, xmax, Mx + 2)(2:Mx+1)';
	 t = linspace(0, T, Mt + 1);
	 dt = T / Mt;


	 ## Auxillary vectors for numerical scheme
	 is = (1:Mx)';
	 r = r0 - r1;
	 A = .5*(sigma^2*is.^2-r*is)*dt;
	 B = -(sigma^2*is.^2 + r0)*dt;
	 C = .5*(sigma^2*is.^2 + r*is)*dt;

	 ## Only Crank - Nicholson so far
	 M1 = spdiags([.5*C, .5*B-1, .5*A], [-1 0 1], Mx, Mx)';
	 M2 = spdiags([-.5*C, -.5*B-1, -.5*A], [-1 0 1], Mx, Mx)';

	 for k = Mt:-1:1
	   V = M2*V;
	   
	   d = [A(1) * left_cond(k) ; zeros(Mx - 2,1) ; C(Mx) * right_cond(k)];
	  
	   V = M1 \ (V + d);
	   
	   V = V .* (x > left_bound(k) & x < right_bound(k));
	   
	 endfor

	 V = V';
endfunction
