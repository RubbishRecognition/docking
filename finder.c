//
//  finder.c
//
//
//  Created by apple on 11.05.16.
//
//
// Uniprot data /Volumes/Basic\ data\ partition\ 1/uniprot_sprot_varsplic.fasta
// BIOGRID data vim /Volumes/Basic\ data\ partition\ 1/BIOGRID-ALL-3.4.136.mitab.txt

//#include "finder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GENE_NAME_LEN 100
#define DEBUG 0

int string_comparing (char* s1, char* s2);

int biogrid_find (char* biogrid_data, char* pattern);

int uniprot_find (char* uniprot_data, char* pattern);


int main (int argc, char* argv[])
{
    //FILE* table;
    
    //able = fopen ("/Volumes/Basic data partition 1/BIOGRID-ALL-3.4.136.mitab.txt", "r");
    
    //char curr_char = 0;
    //fgets(currChar, 2, table);
    //curr_char = fgetc(table);
    /*
     int i = 0;
     for (i = 0; i < 10000; i++)
     {
     printf ("%c", curr_char);
     curr_char = fgetc(table);
     }
     
     fclose (table);
     */

    //printf("kek");
    
    if (argv[1][0] == 'B')
    {
        //work with BIOGRID
        
        biogrid_find (argv[2], argv[3]);
    }
    
    if (argv[1][0] == 'U')
    {
        //work with Uniprot
        
        uniprot_find (argv[2], argv[3]);
    }

    if (argv[1][0] == 'S')
    {
        //work with Uniprot

        string_find (argv[2], argv[3]);
    }
    
    return 0;
}

int string_comparing (char* s1, char* s2)
{
    int flag = 0;
    int len1 = 0;
    int len2 = 0;
    
    while ((flag == 0) & (s1[len1] != 0) & (s2[len2] != 0))
    {
        if (s1[len1] != s2[len2])
            flag = 1;
        else
        {
            len1++;
            len2++;
        }
    }
    
    if ((flag == 0) & (s1[len1] == 0) & (s2[len2] == 0))
        return 1;
    else
        return 0;
}

int biogrid_find (char* biogrid_data, char* pattern)
{
    //printf("kek");
    FILE* data_file = fopen (biogrid_data, "r");
    int ncols = 12;
    char * path = "../tmp/";
    char dest[80] = {};	
    //char dest1[100] = {};
    //char dest2[100] = {};
    //memset (dest, '\0', sizeof(dest));
    //strcpy (dest1, "../tmp/");
    //strcpy (dest2, dest1);
    //FILE* fout = fopen (strncat(path, pattern, 100), "w");
    //FILE* fout = fopen ("../tmp/kek", "w");
    strcat (dest, path);
    strcat (dest, pattern);
    //printf("%s\n", dest); 
    //printf("kek");
    FILE* fout = fopen (dest, "w");
    int curr_col = 0;
    int curr_char = fgetc (data_file);
    
    int is_header = 1;
    
    while (curr_char != EOF)
    {
        curr_char = fgetc (data_file);
        
        if (curr_char == '\t')
        {
            curr_col++;
            curr_char = fgetc (data_file);
            //printf("%c\n", curr_char);
        }
        
        if (curr_char == '\n')
        {
            curr_col = 0;
            is_header = 0;
            curr_char = fgetc (data_file);
            //printf("%c\n", curr_char);
        }
        
        if ((curr_col == 2) & (is_header == 0))
        {
            
            //if (DEBUG)
            //printf("I'm here with sym %c\n", curr_char);
            //column with interactor A
            
            int flag = 0;
            while (flag != 2)
            {
                if (curr_char == ':')
                    flag++;
                
                curr_char = fgetc (data_file);
            }
            
            //now curr_char on first symbol of interactor A
            
            char pattern_A[GENE_NAME_LEN] = {};
            
            //if (DEBUG)
            //printf("Interactor A with sym %c\n", curr_char);
            
            int len = 0;
            while ((curr_char != '\t') & (curr_char != '|'))
            {
                //printf("Interactor A with sym %c\n", curr_char);
                
                pattern_A[len] = curr_char;
                curr_char = fgetc (data_file);
                len++;
            }
            
            while (curr_char != '\t')
                curr_char = fgetc (data_file);
            
            curr_char = fgetc (data_file);
            curr_col++;
            
            //column with interactor B
            
            flag = 0;
            while (flag != 2)
            {
                if (curr_char == ':')
                    flag++;
                
                curr_char = fgetc (data_file);
            }
            
            //now curr_char on first symbol of interactor B
            
            char pattern_B[GENE_NAME_LEN] = {};
            
            len = 0;
            while ((curr_char != '\t') & (curr_char != '|'))
            {
                pattern_B[len] = curr_char;
                curr_char = fgetc (data_file);
                len++;
            }
            
            while (curr_char != '\t')
                curr_char = fgetc (data_file);
            
            curr_char = fgetc (data_file);
            curr_col++;
            
            if (DEBUG)
                printf("Interactor A : %s\n", pattern_A);
            
            if (DEBUG)
                printf("Interactor B : %s\n", pattern_B);
            
            if (string_comparing (pattern, pattern_A))
                fprintf (fout, "%s\n", pattern_B);
            
            if (string_comparing (pattern, pattern_B))
                fprintf (fout, "%s\n", pattern_A);
            
        }
        
    }
    return 0;
}

int uniprot_find (char* uniprot_data, char* pattern)
{
    FILE* data_file = fopen (uniprot_data, "r");
    
    int curr_char = fgetc (data_file);
    
    while (curr_char != EOF)
    {
        int flag = 1;
        while (flag)
        {
            while (curr_char != ' ')
                curr_char = fgetc (data_file);
            curr_char = fgetc (data_file);
            
            if (curr_char == 'G')
            {
                curr_char = fgetc (data_file);
                
                if (curr_char == 'N')
                {
                    curr_char = fgetc (data_file);
                    
                    if (curr_char == '=')
                    {
                        curr_char = fgetc (data_file);
                        flag = 0;
                    }
                }
            }
            
        }
        //now curr_char on first symbol of current Gene
        
        char gene[GENE_NAME_LEN] = {};
        
        int len = 0;
        while (curr_char != '\n')
        {
            gene[len] = curr_char;
            curr_char = fgetc (data_file);
            len++;
        }
        
        curr_char = fgetc (data_file);
        
        if (DEBUG)
            printf("gene : %s\n", gene);
        
        if (string_comparing (pattern, gene))
        {
            while (curr_char != '>')
            {
                printf("%c", curr_char);
                curr_char = fgetc (data_file);
            }
            //printf("\n");
            return 0;
        }
        
        curr_char = fgetc (data_file);
        
    }
    
    return 0;
}

int string_find (char* string_data, char* pattern)
{
    char dest[80] = {};
    FILE* data_file = fopen (string_data, "r");
    //FILE* fout = fopen ("string_output.txt", "w");
    //FILE* fout = fopen (strcpy('string_find_', pattern), "w");
    char * path = "../tmp/";
    //char dest2[100] = {};
    //memset (dest, '\0', sizeof(dest));
    //strcpy (dest, pattern);
    //strcpy (dest2, dest1);
    //FILE* fout = fopen (strncat(path, pattern, 100), "w");
    strcat (dest, path);
    strcat (dest, pattern);
    //printf("%s\n", dest); 
    //printf("kek");
    FILE* fout = fopen (dest, "w");
    //printf("%s\n", strncat(path, pattern, 100));
    //int k = 0;
    //scanf("%d", &k);
    
    int is_header = 1;
    
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    while ((read = getline(&line, &len, data_file)) != -1)
    {
        //printf("here");
        //int k = 0;
        //scanf("%d", &k);
        
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
        
        if (is_header)
        {
            is_header = 0;
            continue;
        }
        else
        {
            char pattern_A[GENE_NAME_LEN] = {};
            char pattern_B[GENE_NAME_LEN] = {};
            
            int n = 0;
            int pnt_flg = 0;
            int gene_len = 0;
            
            while (line[n] != '\t')
            {
                if (pnt_flg)
                {
                    pattern_A[gene_len] = line[n];
                    gene_len += 1;
                }
                
                if (line[n] == '.')
                    pnt_flg = 1;
                
                //n += 1;
                
                //printf("%c, %d, %d, %d", line[n], n, pnt_flg, (line[n] != ' '));
                //int k = 0;
                //scanf("%d", &k);
                
                n += 1;
            }
            
            //printf("%d", n);
            //int k = 0;
            //scanf("%d", &k);
            
            while (line[n] == '\t')
                n += 1;
            
            pnt_flg = 0;
            gene_len = 0;
            
            while (line[n] != '\t')
            {
                if (pnt_flg)
                {
                    pattern_B[gene_len] = line[n];
                    gene_len += 1;
                }
                
                if (line[n] == '.')
                    pnt_flg = 1;
                
                n += 1;
            }
            
            //printf("%d, %s, %s, \n", n, pattern_A, pattern_B);
            //int k = 0;
            //scanf("%d", &k);
            
            if (string_comparing (pattern, pattern_A))
                fprintf (fout, "%s\n", pattern_B);
            
            if (string_comparing (pattern, pattern_B))
                fprintf (fout, "%s\n", pattern_A);
            
        }
        
    }
    
}



















