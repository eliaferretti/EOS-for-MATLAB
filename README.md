# EOS-for-MatLab
Libreria estesa contenente funzioni per il calcolo di grandezze termodinamiche con varie EoS cubiche + Viriale

Prerequisiti EoS 
Dichiarare globali le variabili: 
		-  Tc = Vettore delle temperature critiche [K] 
		-  pc = Vettore delle pressioni critiche [Pa] 
		-  w = Vettore dei fattori acentrici di Pitzer [-] 

Prerequisiti gamma di Wilson 
Dichiarare globale la matrice dei coefficienti di Wilson con il nome “W” (doppia vu maiuscola)
 
Significato dei termini di input 
		-  T/Temp = Temperature [K] 
		-  p/press = Pressione [Pa] 
		-  x = composizione della fase (vettore delle frazioni molari) 
		-  state = Stato fisico (“L” o “V”) 
		-  i/index = indice del composto 

Note sui termini di output 
		-  zeta [-] 
		-  phi [-] 
		-  entalpia residua [J/mol] 
		-  gamma [-] 

Note 
		-  Le funzioni sfruttano le unità di misura del Sistema Internazionale 
		-  Tutte le funzioni sono state testate sui risultati numerici di esercizi proposti nel testo: 
			“Fondamenti di Termodinamica dell’Ingegneria Chimica” – R.Rota 
