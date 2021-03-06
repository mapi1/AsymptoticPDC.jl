"""
    get_sunspot_melanoma_data()

This data is from Andrews and Herzberg. D. F. Andrews, A. M. Herzberg. (1985) Data: A Collection of Problems from Many Fields for the Student and Research Worker. Springer, New York.

Returns `u`:

* `u[:, 1]`: Year 1936-1972 
* `u[:, 2]`: Annual male melanoma incidence (age-adjusted per 10e5) in Connecticut
* `u[:, 3]`: Annual total melanoma incidence (age-adjusted per 10e5) in Connecticut
* `u[:, 4]`: Annual sunspot relative number 
"""
function get_sunspot_melanoma_data()
    u = [1936 1.0 0.9 40
        1937 0.8 0.8 115
        1938 0.8 0.8 100
        1939 1.4 1.3 80
        1940 1.2 1.4 60
        1941 1.0 1.2 40
        1942 1.5 1.7 23
        1943 1.9 1.8 10
        1944 1.5 1.6 10
        1945 1.5 1.5 25
        1946 1.5 1.5 75
        1947 1.6 2.0 145
        1948 1.8 2.5 130
        1949 2.8 2.7 130
        1950 2.5 2.9 80
        1951 2.5 2.5 65
        1952 2.4 3.1 20
        1953 2.1 2.4 10
        1954 1.9 2.2 5
        1955 2.4 2.9 10
        1956 2.4 2.5 60
        1957 2.6 2.6 190
        1958 2.6 3.2 180
        1959 4.4 3.8 175
        1960 4.2 4.2 120
        1961 3.8 3.9 50
        1962 3.4 3.7 35
        1963 3.6 3.3 20
        1964 4.1 3.7 10
        1965 3.7 3.9 15
        1966 4.2 4.1 30
        1967 4.1 3.8 60
        1968 4.1 4.7 105
        1969 4.0 4.4 105
        1970 5.2 4.8 105
        1971 5.3 4.8 80
        1972 5.3 4.8 65]
    return u
end
