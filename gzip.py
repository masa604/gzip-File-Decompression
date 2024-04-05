# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
from huffmantree import HuffmanTree
import numpy as np


class GZIPHeader:
	''' class for reading and storing GZIP header fields '''

	ID1 = ID2 = CM = FLG = XFL = OS = 0
	MTIME = []
	lenMTIME = 4
	mTime = 0

	# bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
	FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0   
	
	# FLG_FTEXT --> ignored (usually 0)
	# if FLG_FEXTRA == 1
	XLEN, extraField = [], []
	lenXLEN = 2
	
	# if FLG_FNAME == 1
	fName = ''  # ends when a byte with value 0 is read
	
	# if FLG_FCOMMENT == 1
	fComment = ''   # ends when a byte with value 0 is read
		
	# if FLG_HCRC == 1
	HCRC = []
		
		
	
	def read(self, f):
		''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

		# ID 1 and 2: fixed values
		self.ID1 = f.read(1)[0]  
		if self.ID1 != 0x1f: return -1 # error in the header
			
		self.ID2 = f.read(1)[0]
		if self.ID2 != 0x8b: return -1 # error in the header
		
		# CM - Compression Method: must be the value 8 for deflate
		self.CM = f.read(1)[0]
		if self.CM != 0x08: return -1 # error in the header
					
		# Flags
		self.FLG = f.read(1)[0]
		
		# MTIME
		self.MTIME = [0]*self.lenMTIME
		self.mTime = 0
		for i in range(self.lenMTIME):
			self.MTIME[i] = f.read(1)[0]
			self.mTime += self.MTIME[i] << (8 * i) 				
						
		# XFL (not processed...)
		self.XFL = f.read(1)[0]
		
		# OS (not processed...)
		self.OS = f.read(1)[0]
		
		# --- Check Flags
		self.FLG_FTEXT = self.FLG & 0x01
		self.FLG_FHCRC = (self.FLG & 0x02) >> 1
		self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
		self.FLG_FNAME = (self.FLG & 0x08) >> 3
		self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4
					
		# FLG_EXTRA
		if self.FLG_FEXTRA == 1:
			# read 2 bytes XLEN + XLEN bytes de extra field
			# 1st byte: LSB, 2nd: MSB
			self.XLEN = [0]*self.lenXLEN
			self.XLEN[0] = f.read(1)[0]
			self.XLEN[1] = f.read(1)[0]
			self.xlen = self.XLEN[1] << 8 + self.XLEN[0]
			
			# read extraField and ignore its values
			self.extraField = f.read(self.xlen)
		
		def read_str_until_0(f):
			s = ''
			while True:
				c = f.read(1)[0]
				if c == 0: 
					return s
				s += chr(c)
		
		# FLG_FNAME
		if self.FLG_FNAME == 1:
			self.fName = read_str_until_0(f)
		
		# FLG_FCOMMENT
		if self.FLG_FCOMMENT == 1:
			self.fComment = read_str_until_0(f)
		
		# FLG_FHCRC (not processed...)
		if self.FLG_FHCRC == 1:
			self.HCRC = f.read(2)
			
		return 0
			



class GZIP:
    ''' class for GZIP decompressing file (if compressed with deflate) '''

    gzh = None
    gzFile = ''
    fileSize = origFileSize = -1
    numBlocks = 0
    f = None

    bits_buffer = 0
    available_bits = 0

    def __init__(self, filename):
        self.gzFile = filename
        self.f = open(filename, 'rb')
        self.f.seek(0,2)
        self.fileSize = self.f.tell()
        self.f.seek(0)

    def decompress(self):
        ''' main function for decompressing the gzip file with deflate algorithm '''

        numBlocks = 0

        # get original file size: size of file before compression
        origFileSize = self.getOrigFileSize()
        print(origFileSize)

        # read GZIP header
        error = self.getHeader()
        if error != 0:
            print('Formato invalido!')
            return

        # show filename read from GZIP header
        print(self.gzh.fName)

        # MAIN LOOP - decode block by block
        BFINAL = 0    
        while not BFINAL == 1:    

            BFINAL = self.readBits(1)

            BTYPE = self.readBits(2)                  
            if BTYPE != 2:
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks+1))
                return

            # --- STUDENTS --- ADD CODE HERE
            # 
            # ex 1
            HLIT = 257+self.readBits(5)
            HDIST = 1+self.readBits(5)
            HCLEN = 4+self.readBits(4)
            print("HLIT =",HLIT)
            print("HDIST =",HDIST)
            print("HCLEN =",HCLEN)
            
            # ex 2
            vetorOrdem = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
            n_array=self.ex2(HCLEN, vetorOrdem)
            print(n_array)
            
            # ex 3
            codigosArvore=self.ex3(n_array)
            print(codigosArvore)
            tree=self.transformaArrTree(codigosArvore)
            
            # ex 4 
            HLITarr=self.ex4e5(tree, HLIT)
            #print("HLITarr: ",list(HLITarr))
            
            # ex 5
            HDISTarr=self.ex4e5(tree, HDIST)
            #print("HDISTarr: ",list(HDISTarr))
            
            # ex 6
            arrLIT,arrDIST,treeLIT,treeDIST=self.ex6(HLITarr,HDISTarr)
            #print("arrLIT: ",arrLIT)
            #print("arrDIST: ",arrDIST)
            
            # ex7
            vetor,texto=self.exercicio7(treeLIT, treeDIST)
            print(texto)
            
            #ex 8
            flname=GZIPHeader.fName
            flname=str(flname)
            with open(flname,'w') as ficheiro:
                print("Ficheiro",flname,"criado com sucesso")
                ficheiro.write(texto)

            # update number of blocks read
            numBlocks += 1

        # close file            
        self.f.close()    
        print("End: %d block(s) analyzed." % numBlocks)
    
    
    def getOrigFileSize(self):
        ''' reads file size of the original file (before compression) - ISIZE '''
        
        # saves the current position of the file pointer
        fp = self.f.tell()
        
        # jumps to the end-4 position
        self.f.seek(self.fileSize-4)
        
        # reads the last 4 bytes (LITTLE ENDIAN)
        sz = 0
        for i in range(4): 
            sz += self.f.read(1)[0] << (8*i)
        
        # restores the file pointer to its original position
        self.f.seek(fp)
        
        return sz
    
    def getHeader(self):  
        ''' reads GZIP header'''

        self.gzh = GZIPHeader()
        header_error = self.gzh.read(self.f)
        return header_error

    def readBits(self, n, keep=False):
        ''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

        while n > self.available_bits:
            self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
            self.available_bits += 8
        
        mask = (2**n)-1
        value = self.bits_buffer & mask

        if not keep:
            self.bits_buffer >>= n
            self.available_bits -= n

        return value
    
    def ex2(self,HCLEN,vetorOrdem):
        n_array=np.zeros(len(vetorOrdem),dtype=np.uint16)
        for i in range(HCLEN):
            n_array[vetorOrdem[i]]=self.readBits(3)
        return n_array
    
    def ex3(self,n_array):
        # Passo 1 (doc2) Conte o número de códigos para cada comprimento de código. Seja bl_count[N] o número de códigos de comprimento N, N >= 1.
        numeros = np.unique(n_array)
        numeros = np.sort(numeros)
        print(numeros)  # Coloca no array os elementos uma só vez (evita repetições)
        aux = max(numeros)

        bl_count = np.bincount(n_array)  # Conta o número de ocorrências de cada número em n_array
        print(bl_count)

        # Passo 2
        code = 0
        bl_count[0] = 0
        next_code = np.zeros(aux + 1, int)

        for bits in range(1, aux + 1):
            code = (code + bl_count[bits - 1]) << 1
            next_code[bits] = code

        print(next_code)
        # Passo 3
        tree = []
        for i in range(len(n_array)):
            tree.append("")  # Torna os elementos como strings

        for n in range(len(n_array)):
            len_value = n_array[n]  # len_value é atribuído ao valor da lista n_array no índice atual
            if len_value != 0:
                tree[n] = np.binary_repr(next_code[len_value], width=n_array[
                    n])  # representação binária do número de entrada como uma string, do tamanho de n_array[n]
                next_code[len_value] += 1
        
        return tree
    
    def transformaArrTree(self,arr):
        HuffmanCode = HuffmanTree()
        for i in range (len(arr)) :
            if arr[i] != '':
                HuffmanCode.addNode(arr[i],i,True)
        return HuffmanCode
    
    def ex4e5(self,hufftree,tamanho):
        hufftree.resetCurNode()
        arr=[0 for i in range(tamanho)]
        i=0
        while i<tamanho:
            no=self.lerArvore(hufftree)
            if no==16:
                for j in range(3+self.readBits(2)):
                    arr[i]=arr[i-1]
                    i+=1
            elif no == 17:
                for j in range(3 + self.readBits(3)):
                    arr[i] = 0
                    i += 1
            elif no == 18:
                repeat = 11 + self.readBits(7)
                for j in range(repeat):
                    arr[i] = 0
                    i += 1
            else:
                arr[i] = no
                i += 1
        
        return arr
            
    def lerArvore(self,tree):
        tree.resetCurNode()
        while True:
            proxBit=self.readBits(1)
            no=tree.nextNode(str(proxBit))
            if no==-1:
                print("Arvore de huffman mal construida!")
                exit()
            elif no!=-2:
                return no
            
    def ex6(self,literais,distancias):
        arrLit=self.ex3(literais)                   # cria um array com os códigos de Huffman dos literais
        arrDist=self.ex3(distancias)                # cria um array com os códigos Huffman das distancias/comprimentos
        TreeLIT=self.transformaArrTree(arrLit)      # transforma o array numa arvore com os codigos dos literais
        TreeDIST=self.transformaArrTree(arrDist)    # transforma o array numa arvore com os codigos das distancias
        
        return arrLit,arrDist,TreeLIT,TreeDIST
    
    def exercicio7(self,treeLIT,treeDIST):
        strLIT = ""
        strDIST = ""
        texto = ""
        vetor = []
        codeLIT = treeLIT.findNode(strLIT)
        codeDIST = treeDIST.findNode(strDIST)
        
        
        while (codeLIT != 256):
            while (codeLIT == -2):
                strLIT += bin(self.readBits(1))[2:]
                codeLIT = treeLIT.findNode(strLIT)
                
            if (codeLIT < 256):
                vetor.append(codeLIT)
                texto += chr(codeLIT)
                
            
            elif (codeLIT == 256):
                break
            else:
                if (257 <= codeLIT <= 264):
                    comp = (codeLIT - 257) + 3
                elif (265 <= codeLIT <= 268):
                    comp = 2*(codeLIT - 265) + 11 + self.readBits(1)
                elif (269 <= codeLIT <= 272):
                    comp = 4*(codeLIT - 269) + 19 + self.readBits(2)
                elif (273 <= codeLIT <= 276):
                    comp = 8*(codeLIT - 273) + 35 + self.readBits(3)
                elif 277 <= codeLIT <= 280:
                    comp = 16*(codeLIT - 277) + 67 + self.readBits(4)
                elif 281 <= codeLIT <= 284:
                    comp = 32*(codeLIT - 281) + 115 + self.readBits(5)
            
                        
                while (codeDIST == -2):
                    strDIST += bin(self.readBits(1))[2:]
                    codeDIST = treeDIST.findNode(strDIST)
                    
                if (codeDIST <= 3):    
                    recuar = codeDIST + 1
                elif(codeDIST == 4):
                     recuar = 5 + self.readBits(1)
                elif(codeDIST == 5):
                    recuar = 7 + self.readBits(1)
                elif(codeDIST == 6):
                    recuar = 9 + self.readBits(2)
                elif(codeDIST == 7):
                    recuar = 13 + self.readBits(2) 
                elif(codeDIST == 8):
                    recuar = 17 + self.readBits(3)
                elif(codeDIST == 9):
                    recuar = 25 + self.readBits(3)
                elif(codeDIST == 10):
                    recuar = 33 + self.readBits(4)
                elif(codeDIST == 11):
                    recuar = 49 + self.readBits(4)
                elif(codeDIST == 12):
                    recuar = 65 + self.readBits(5)
                elif(codeDIST == 13):
                    recuar = 97 + self.readBits(5)
                elif(codeDIST == 14):
                    recuar = 129 + self.readBits(6)
                elif(codeDIST == 15):
                    recuar = 193 + self.readBits(6)
                elif(codeDIST == 16):
                    recuar = 257 + self.readBits(7)
                elif(codeDIST == 17):
                    recuar = 385 + self.readBits(7)
                elif(codeDIST == 18):
                    recuar = 513 + self.readBits(8)
                elif(codeDIST == 19):
                    recuar = 769 + self.readBits(8)
                elif(codeDIST == 20):
                    recuar = 1025 + self.readBits(9)
                elif(codeDIST == 21):
                    recuar = 1537 + self.readBits(9)
                elif(codeDIST == 22):
                    recuar = 2049 + self.readBits(10)
                elif(codeDIST == 23):
                    recuar = 3073 + self.readBits(10)
                elif(codeDIST == 24):
                    recuar = 4097 + self.readBits(11)
                elif(codeDIST == 25):
                    recuar = 6145 + self.readBits(11)
                elif(codeDIST == 26):
                    recuar = 8193 + self.readBits(12)
                elif(codeDIST == 27):
                    recuar = 12289 + self.readBits(12)
                elif(codeDIST == 28):
                    recuar = 16385 + self.readBits(13)
                elif(codeDIST == 29):
                    recuar = 24577 + self.readBits(13)
        
                for i in range(comp):
                    vetor.append(vetor[-int(recuar)])
                    texto += chr(vetor[-1])
                            
            strLIT = ""
            strDIST = ""
                
            codeLIT = treeLIT.findNode(strLIT)
            codeDIST = treeDIST.findNode(strDIST)
         
        return vetor, texto
        

if __name__ == '__main__':
    # gets filename from command line if provided
    fileName="FAQ.txt.gz"
    if len(sys.argv)>1:
        fileName=sys.argv[1]	

    GZIPHeader.fName="FAQ.txt"
    gz=GZIP(fileName)
    gz.decompress()		
	
