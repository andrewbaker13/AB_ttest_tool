// Constants and utility functions
const CREATED_DATE = new Date('2025-11-06').toLocaleDateString();
let modifiedDate = new Date().toLocaleDateString();

// Error function approximation (Abramowitz and Stegun)
function erf(x) {
    const sign = Math.sign(x);
    x = Math.abs(x);
    
    const a1 = 0.254829592;
    const a2 = -0.284496736;
    const a3 = 1.421413741;
    const a4 = -1.453152027;
    const a5 = 1.061405429;
    const p = 0.3275911;

    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    
    return sign * y;
}

function normCdf(x) {
    return 0.5 * (1 + erf(x / Math.sqrt(2)));
}

function tCdf(t, df) {
    // Approximation of Student's t CDF using normal CDF
    // More accurate for df > 30, but reasonable for visualization
    return normCdf(t);
}

// Calculate effect size (Cohen's d)
function calculateCohensD(mean1, sd1, n1, mean2, sd2, n2) {
    // Pooled standard deviation
    const pooledSD = Math.sqrt(((n1 - 1) * sd1 * sd1 + (n2 - 1) * sd2 * sd2) / (n1 + n2 - 2));
    return Math.abs(mean1 - mean2) / pooledSD;
}

// Calculate statistical power
function calculatePower(d, n1, n2, alpha = 0.05) {
    const ncp = d * Math.sqrt((n1 * n2) / (n1 + n2)); // Non-centrality parameter
    const criticalT = -normCdf(alpha/2);
    // Approximation of power using normal distribution
    return 1 - normCdf(criticalT - ncp) + normCdf(-criticalT - ncp);
}

function getEffectSizeInterpretation(d) {
    if (d < 0.2) return "Very small effect size";
    if (d < 0.5) return "Small effect size";
    if (d < 0.8) return "Medium effect size";
    return "Large effect size";
}

// Calculate Welch's t-test statistics
function calculateWelchTTest(mean1, sd1, n1, mean2, sd2, n2, delta0 = 0) {
    const diff = mean1 - mean2 - delta0;
    const se1 = (sd1 * sd1) / n1;
    const se2 = (sd2 * sd2) / n2;
    const se = Math.sqrt(se1 + se2);
    
    const t = diff / se;
    
    // Welch–Satterthwaite degrees of freedom
    const df = Math.pow(se1 + se2, 2) / 
        (Math.pow(se1, 2) / (n1 - 1) + Math.pow(se2, 2) / (n2 - 1));
    
    // Calculate effect size and power
    const cohensD = calculateCohensD(mean1, sd1, n1, mean2, sd2, n2);
    const power = calculatePower(cohensD, n1, n2);
    
    return { t, df, se, cohensD, power };
}

// Update visualization and interpretation
function updateResults() {
    const mean1 = parseFloat(document.getElementById('mean1').value);
    const sd1 = parseFloat(document.getElementById('sd1').value);
    const n1 = parseInt(document.getElementById('n1').value);
    const mean2 = parseFloat(document.getElementById('mean2').value);
    const sd2 = parseFloat(document.getElementById('sd2').value);
    const n2 = parseInt(document.getElementById('n2').value);
    const delta0 = parseFloat(document.getElementById('delta0').value) || 0;
    const alpha = parseFloat(document.getElementById('alpha').value);

    // Validate inputs
    if ([mean1, sd1, n1, mean2, sd2, n2].some(isNaN) || sd1 <= 0 || sd2 <= 0 || n1 < 2 || n2 < 2) {
        document.getElementById('interpretation').textContent = 
            'Please enter valid values. Standard deviations must be positive and sample sizes must be at least 2.';
        return;
    }

    // Calculate test statistics
    const { t, df, se, cohensD, power } = calculateWelchTTest(mean1, sd1, n1, mean2, sd2, n2, delta0);
    const criticalT = -normCdf(alpha/2); // Approximation
    const pValue = 2 * (1 - normCdf(Math.abs(t)));
    
    // Calculate confidence interval
    const ciMargin = criticalT * se;
    const ciLower = (mean1 - mean2) - ciMargin;
    const ciUpper = (mean1 - mean2) + ciMargin;

    // Create visualization
    createPlot(mean1, sd1, n1, mean2, sd2, n2, delta0, criticalT);

    // Update interpretation
    updateInterpretation(t, df, pValue, ciLower, ciUpper, delta0, alpha, cohensD, power);

    // Update modified date
    modifiedDate = new Date().toLocaleDateString();
    document.getElementById('modified-date').textContent = modifiedDate;
}

function createPlot(mean1, sd1, n1, mean2, sd2, n2, delta0, criticalT) {
    // Calculate confidence intervals
    const alpha = parseFloat(document.getElementById('alpha').value);
    const ci1Margin = criticalT * (sd1 / Math.sqrt(n1));
    const ci2Margin = criticalT * (sd2 / Math.sqrt(n2));
    
    // Calculate overall range for x-axis
    const minX = Math.min(mean1 - 3*sd1/Math.sqrt(n1), mean2 - 3*sd2/Math.sqrt(n2));
    const maxX = Math.max(mean1 + 3*sd1/Math.sqrt(n1), mean2 + 3*sd2/Math.sqrt(n2));
    
    // Group Comparison Plot
    const trace1 = {
        x: [mean1],
        y: ['Group 1'],
        error_x: {
            type: 'data',
            array: [ci1Margin],
            visible: true,
            color: 'rgb(31, 119, 180)'
        },
        type: 'scatter',
        mode: 'markers',
        marker: {
            color: 'rgb(31, 119, 180)',
            size: 10
        },
        name: `Group 1 (n=${n1})`
    };

    const trace2 = {
        x: [mean2],
        y: ['Group 2'],
        error_x: {
            type: 'data',
            array: [ci2Margin],
            visible: true,
            color: 'rgb(255, 127, 14)'
        },
        type: 'scatter',
        mode: 'markers',
        marker: {
            color: 'rgb(255, 127, 14)',
            size: 10
        },
        name: `Group 2 (n=${n2})`
    };

    // Difference Plot
    const diffMean = mean1 - mean2;
    const diffSE = Math.sqrt((sd1*sd1/n1) + (sd2*sd2/n2));
    const diffCI = criticalT * diffSE;

    const traceDiff = {
        x: [diffMean],
        y: ['Difference'],
        error_x: {
            type: 'data',
            array: [diffCI],
            visible: true,
            color: 'rgb(44, 160, 44)'
        },
        type: 'scatter',
        mode: 'markers',
        marker: {
            color: 'rgb(44, 160, 44)',
            size: 10
        },
        name: 'Difference'
    };

    const nullLine = {
        x: [delta0, delta0],
        y: ['Difference', 'Difference'],
        type: 'scatter',
        mode: 'lines',
        line: {
            color: 'red',
            dash: 'dash'
        },
        name: 'Null Hypothesis'
    };

    // Create subplot layout
    const layout = {
        grid: {
            rows: 2,
            columns: 1,
            pattern: 'independent',
            roworder: 'top to bottom'
        },
        height: 600,
        title: 'Group Comparison and Difference Analysis',
        showlegend: true,
        legend: { orientation: 'h', y: -0.2 },
        annotations: [],
        xaxis: {
            title: 'Value',
            range: [minX, maxX]
        },
        xaxis2: {
            title: 'Difference (Group 1 - Group 2)',
            range: [diffMean - 3*diffSE, diffMean + 3*diffSE]
        }
    };

    // Add interpretation annotations
    const pValue = 2 * (1 - normCdf(Math.abs((diffMean - delta0)/diffSE)));
    const significant = pValue < alpha;
    
    layout.annotations.push({
        x: mean1,
        y: 'Group 1',
        xref: 'x',
        yref: 'y',
        text: `Mean: ${mean1.toFixed(2)}`,
        showarrow: false,
        yshift: 20
    }, {
        x: mean2,
        y: 'Group 2',
        xref: 'x',
        yref: 'y',
        text: `Mean: ${mean2.toFixed(2)}`,
        showarrow: false,
        yshift: 20
    }, {
        x: diffMean,
        y: 'Difference',
        xref: 'x2',
        yref: 'y2',
        text: `Diff: ${diffMean.toFixed(2)}\n${significant ? 'Significant' : 'Not Significant'}`,
        showarrow: false,
        yshift: 20,
        font: {
            color: significant ? 'rgb(44, 160, 44)' : 'rgb(214, 39, 40)'
        }
    });

    Plotly.newPlot('visualization', 
        // Split traces between subplots using different domains
        [trace1, trace2, traceDiff, nullLine], 
        layout
    );
}

function updateInterpretation(t, df, pValue, ciLower, ciUpper, delta0, alpha, cohensD, power) {
    const interpretation = document.getElementById('interpretation');
    
    const roundedT = t.toFixed(3);
    const roundedDf = df.toFixed(1);
    const roundedP = pValue.toFixed(4);
    const roundedCiLower = ciLower.toFixed(3);
    const roundedCiUpper = ciUpper.toFixed(3);
    const roundedD = cohensD.toFixed(3);
    const roundedPower = (power * 100).toFixed(1);
    
    let text = `
        <h3>Statistical Results:</h3>
        <div class="results-grid">
            <div class="result-item">
                <h4>Primary Test Statistics</h4>
                <p><strong>Test Statistic:</strong> t = ${roundedT} (df = ${roundedDf})</p>
                <p><strong>P-value:</strong> ${roundedP}</p>
                <p><strong>${(1-alpha)*100}% Confidence Interval:</strong> (${roundedCiLower}, ${roundedCiUpper})</p>
            </div>
            
            <div class="result-item">
                <h4>Effect Size Analysis</h4>
                <p><strong>Cohen's d:</strong> ${roundedD}</p>
                <p><strong>Interpretation:</strong> ${getEffectSizeInterpretation(cohensD)}</p>
                <p><strong>Statistical Power:</strong> ${roundedPower}%</p>
            </div>
        </div>
        
        <h3>Detailed Interpretation:</h3>
        <div class="interpretation-details">
            <h4>Statistical Significance:</h4>
            <p>${
                pValue < alpha 
                    ? `We reject the null hypothesis that the true difference equals ${delta0}.`
                    : `We fail to reject the null hypothesis that the true difference equals ${delta0}.`
            }</p>
            <p>The data ${
                pValue < alpha 
                    ? 'provides' 
                    : 'does not provide'
            } sufficient evidence of a difference from the hypothesized value at the ${alpha} significance level.</p>
            
            <h4>Practical Significance:</h4>
            <p>The effect size (Cohen's d = ${roundedD}) indicates ${getEffectSizeInterpretation(cohensD).toLowerCase()}. 
               This means the difference between the groups is ${
                   cohensD < 0.2 ? 'very small and might not be practically meaningful' :
                   cohensD < 0.5 ? 'small but might be meaningful in some contexts' :
                   cohensD < 0.8 ? 'moderate and likely practically meaningful' :
                   'large and practically significant'
               }.</p>
            
            <h4>Statistical Power:</h4>
            <p>The test has ${roundedPower}% power to detect the observed effect size. ${
                power < 0.8 
                    ? 'This is below the conventional 80% threshold, suggesting the test might be underpowered. Consider increasing sample sizes.'
                    : 'This exceeds the conventional 80% threshold, indicating adequate power to detect the observed effect.'
            }</p>
            
            <div class="educational-note">
                <h4>📚 Learning Note:</h4>
                <p>Remember that statistical significance (p-value) and practical significance (effect size) tell different stories:
                <ul>
                    <li>P-value tells us how likely we would observe such results under the null hypothesis</li>
                    <li>Effect size tells us about the magnitude of the difference, regardless of sample size</li>
                    <li>Power tells us about our ability to detect true effects when they exist</li>
                </ul>
                </p>
            </div>
        </div>
    `;
    
    interpretation.innerHTML = text;
}

// Event listeners
document.addEventListener('DOMContentLoaded', () => {
    // Set dates
    document.getElementById('created-date').textContent = CREATED_DATE;
    document.getElementById('modified-date').textContent = modifiedDate;
    
    // Add input listeners
    const inputs = [
        'mean1', 'sd1', 'n1',
        'mean2', 'sd2', 'n2',
        'delta0', 'alpha'
    ];
    
    inputs.forEach(id => {
        document.getElementById(id).addEventListener('input', updateResults);
    });
    
    // Initial update
    updateResults();
});