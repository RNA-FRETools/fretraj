# Making FRET predictions

Once a donor and acceptor ACV is calculated, we can predict a **mean transfer efficiency** $\langle E\rangle$ based on the distance between the two volumes.

$$ \langle\ E\rangle = \frac{1}{n\,m}\sum_{i=1}^n\sum_{j=1}^m\frac{1}{1+\parallel\mathbf{R}_{A,j}-\parallel\mathbf{R}_{D,i}\parallel^6 / R_0^6} $$

Here, $\mathbf{R}_{D,i}$ and $\mathbf{R}_{A,j}$ are the coordinate vectors of the points $i$ and $j$ in the donor and acceptor volume {cite}`Sindbert.2011, Hellenkamp.2018`.

1. **R0**: Förster radius of the dye pair
2. **# distances**: number of distances to sample for the FRET calculation. The algorithm chooses *n* pairs of randomly distributed points in the donor and acceptor cloud and calculates their distance
3. choose the **donor** and **acceptor** cloud from the drop down
4. **Calculate FRET**: start the FRET calculation. The mean FRET efficiency $\langle E\rangle$, the mean inter-dye distance $\langle R_{DA}\rangle$ as well as the distances between the mean dye positions $\langle R_{MP}\rangle$ or the two attachment sites are displayed in the table. The FRET parameters are also saved to a JSON file.

```{note}  
The style of the ACV clouds can be tuned with the *ACV visualization* settings:
- **contour levels** of the accessible and contact volume
- **b-factor** and the **gaussian resolution** define the smoothness of the cloud
- **grid buffer** around the ACV
- **transparency** of the cloud
```
