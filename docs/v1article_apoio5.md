

Com base nos documentos fornecidos, aqui está uma introdução estruturada para um artigo científico de alto impacto, abordando a avaliação espectral, as limitações das técnicas de observação e as soluções emergentes, acompanhada das referências bibliográficas correspondentes.

---

**Introdução**

As ondas de superfície oceânica desempenham um papel fundamental na interface oceano-atmosfera, sendo responsáveis por mediar as trocas de energia, momento e massa entre os dois meios. Historicamente, a descrição do estado de mar baseou-se em parâmetros integrais ( *bulk* ), como a altura significativa ($H_s$) e o período médio ($T_m$). No entanto, a maioria das regiões oceânicas é caracterizada por condições multimodais, onde sistemas de mar de vento coexistem com múltiplos trens de ondulação ( *swells* ) provenientes de fontes distintas. Nestes cenários complexos, os parâmetros integrais tornam-se fisicamente ambíguos, falhando em capturar a estrutura energética detalhada necessária para aplicações que vão da engenharia oceânica à modelagem climática global.

Neste contexto, o Radar de Abertura Sintética (SAR) consolidou-se como o único sensor capaz de fornecer informações bidimensionais sobre a superfície oceânica em escala global e com alta resolução espacial, independentemente da cobertura de nuvens. No entanto, a recuperação do espectro de ondas a partir de imagens SAR é inerentemente complexa devido ao mapeamento não-linear entre o oceano e a imagem. O principal mecanismo de distorção é o  *velocity bunching* , causado pelo movimento orbital das ondas curtas, que introduz um corte espectral na direção do voo do satélite ( *azimuth cut-off* ). Este efeito atua como um filtro passa-baixa, eliminando informações de comprimentos de onda inferiores a 200 m e limitando a observação direta ao campo de  *swell* . Adicionalmente, espectros derivados de imagens de intensidade sofrem de uma ambiguidade inerente de 180 graus na direção de propagação.

Para superar a insuficiência dos parâmetros totais e avaliar a precisão de modelos numéricos de terceira geração, como o WaveWatch III (WW3), a técnica de particionamento espectral tornou-se indispensável. O particionamento permite decompor o espectro em sistemas de ondas individuais, cada um associado a um evento meteorológico específico. Contudo, a aplicação prática revela gargalos significativos: discrepâncias sistemáticas entre as partições observadas pelo SAR e as previstas pelos modelos, frequentemente agravadas pela geração de partições espúrias ou artificiais devido ao ruído instrumental. Além disso, a atribuição cruzada ( *cross-assignment* ) de partições entre observação e modelo permanece um desafio técnico em mares cruzados.

Soluções contemporâneas têm focado em metodologias híbridas e no uso de novos descritores espectrais. O emprego de espectros cruzados ( *cross-spectra* ) permitiu mitigar a ambiguidade direcional e reduzir o ruído de  *speckle* . Recentemente, o parâmetro *Mean rAnge Cross-Spectrum* (MACS) emergiu como uma ferramenta robusta para validar assinaturas espectrais sem a necessidade de esquemas complexos de inversão total. Paralelamente, o uso de algoritmos de aprendizado de máquina, como Redes Neurais e *eXtreme Gradient Boosting* (XGBoost), tem demonstrado eficácia na calibração de parâmetros particionados e na criação de sinalizadores de qualidade ( *Quality Flags* ) automatizados. Esses avanços possibilitam a filtragem de dados corrompidos por não-linearidades, provendo uma base de dados mais precisa para o monitoramento contínuo do clima de ondas global.

---

### Referências Bibliográficas

* **Ardhuin, F., et al. (2019).** Observing Sea States.  *Frontiers in Marine Science* , 6:124.
* **Benchaabane, A., et al. (2025).** Developing a Quality Flag for SAR Ocean Wave Spectrum Partitioning with Machine Learning.  *Remote Sensing* , 17, 3191.
* **Chapron, B., Johnsen, H., & Garello, R. (2001).** Wave and wind retrieval from SAR images of the ocean.  *Annals of Telecommunications* , 56, 682–699.
* **Engen, G., & Johnsen, H. (1995).** SAR-Ocean Wave Inversion Using Image Cross Spectra.  *IEEE Transactions on Geoscience and Remote Sensing* , 33(4), 1047-1056.
* **Hasselmann, K., & Hasselmann, S. (1991).** On the nonlinear mapping of an ocean wave spectrum into a synthetic aperture radar image spectrum and its inversion.  *Journal of Geophysical Research* , 96(C6), 10713–10729.
* **Hasselmann, K., et al. (1985).** Theory of Synthetic Aperture Radar Ocean Imaging: A MARSEN View.  *Journal of Geophysical Research* , 90(C3), 4659–4686.
* **Hasselmann, S., et al. (1996).** An improved algorithm for the retrieval of ocean wave spectra from synthetic aperture radar image spectra.  *Journal of Geophysical Research* , 101(C7), 16615–16629.
* **Li, H., et al. (2021).** Global ocean SAR cross-spectral parameter MACS.  *Journal of Geophysical Research: Oceans* .
* **Li, X.-M., Lehner, S., & Bruns, T. (2011).** Ocean wave integral parameter measurements using Envisat ASAR wave mode data.  *IEEE Transactions on Geoscience and Remote Sensing* , 49(1), 155–174.
* **Portilla-Yandún, J., & Bidlot, J.-R. (2025).** A Global Ocean Spectral Wave Climate Based on ERA-5 Data: GLOSWAC-5.  *Journal of Geophysical Research: Oceans* , 130, e2025JC022629.
* **Stopa, J. E., & Mouche, A. (2017).** Significant wave heights from Sentinel-1 SAR: Validation and applications.  *Journal of Geophysical Research: Oceans* , 122, 1827–1848.
* **Wang, X., Wang, X., & Ge, L. (2023).** Validation and calibration of partitioned integral ocean wave parameters from co-polarized synthetic aperture radar data.  *Remote Sensing of Environment* , 287, 113463.
