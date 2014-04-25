unsigned long long ComputeBinomialCoefficient( int n, int k )
{
        // Run-time assert to ensure correct behavior
  //     assert( n > k && n > 1 ); esta linea es para gente muy cuidadosa

        // Exploit the symmetry in the line x = k/2:
        if( k > n - k )
                k = n - k;

        unsigned long long c(1);
        // Perform the product over the space i = [1...k]
        for( int i = 1; i < k+1; i++ )
        {
                c *= n - (k - i);
                c /= i;
        }

        return c;

};


