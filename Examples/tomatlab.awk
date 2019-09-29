{ 
   if (NR==1) {
	printf("J=zeros(%d,%d);\n",$1,$2); 
   }
   else if ($1!=0) {
	printf("J(%d,%d)=%s;\n",$1,$2,$3);
   }
}
END { 
	printf("J=sparse(J);\n");
    }
