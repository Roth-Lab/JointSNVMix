# cython: profile=True
import string

cdef extern from "string.h":
    char * strtok( char * , char * )
    
cdef extern from "string.h":
    char * strncpy( char * , char * , size_t )
	
cdef extern from "stdlib.h":
    void free( void * ptr )
    void * malloc( size_t size )
    int atoi ( char * )
    
cdef extern from "ctype.h":
    int toupper( int )

cdef extern from "ctype.h":
    int isdigit( int )

cpdef list parse_call_string( char * ref_base, char * call_string ):
    cdef int i, j, digit_length

    cdef int skip_value = 0

    cdef list bases = []
    
    cdef char call_char
    
    cdef char * digit
    
    i = -1

    for call_char in call_string:
        call_char = toupper( < int > call_char )
        
        i += 1   
        
        if skip_value > 0:
            skip_value -= 1
            continue
	
        if call_char in ( b',', b'.' ):
            bases.append( ref_base )

        elif call_char in ( b'A', b'C', b'T', b'G' ):
            bases.append( < bytes > call_char )

        elif call_char in ( b"$", b'N', b'*' ):
            pass

        # End of read skip it and next value which holds no call information.
        elif call_char == b'^':
            skip_value = 1

        # Insertion/deletion info. Skip this the next one and all character relating.
        elif call_char in ( b'+', b'-' ):
            j = i + 1                

            while isdigit( call_string[j] ):
                j += 1

            skip_value += j - i - 1
            
            digit_length = ( j - i - 1 )
            
            digit = < char *> malloc( ( digit_length + 1 ) * sizeof( char ) )
            
            strncpy( digit, call_string + ( i + 1 ), digit_length )
            digit[digit_length] = "\0"
            
            skip_value += atoi( digit )
            
            free( digit )
        else:
            print "This should not happen"
            raise Exception( 'Unparasable char {0}'.format( call_char ) )

    return bases
